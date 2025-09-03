function mu = score_pose_contrast_whole_cylinder_shell(C,u,rad_mm,len_mm,Avox2ras,ctVol,voxdims, insulating_spacer_length)
% cylinder frame
% \text{point}(t,r,\theta)=C+t\,u+r(\cos\theta\,v_1+\sin\theta\,v_2)
% improvement suggestion: less rely on predefined rad_mm and len_mm, 
% in the out model, we fitted the bolb and get rmm_975, which is the rad_mm
% that encompass 97.5% voxels, we will pass it to rad_mm
    % ---- input hygiene ----
    C = reshape(double(C),1,3);
    u = reshape(double(u),1,3);
    u = u./max(norm(u),eps);

    [v1,v2] = orthobasis_row(u);   % 1x3 each

    % Sampling density tied to native voxel sizes
    step_mm = max(min(voxdims)/2, 0.2);
    t  = -len_mm/5*3 : step_mm : +len_mm/5*3; % introduce slightly overlap
    t_shell = (-(len_mm/2 + insulating_spacer_length)) : step_mm : (+(len_mm/2 + insulating_spacer_length));


    % check if C at and u are correct - yes
    do_viz = false;              % <— flip to false to disable
    thrHU  = 3000;              % threshold for “metal”
    vizR   = 60;                % half-width of crop in mm
    vecL   = 10;                % arrow length in mm
    
    if do_viz
        [m,n,p] = size(ctVol);
        [Xm, Ym, Zm] = meshgrid(1:n, 1:m, 1:p);
        figure('Color','w'); hold on;
        BW = ctVol > thrHU;
        h_metal = patch(isosurface(Xm, Ym, Zm, BW, 0.5));
        % make 3d surface look smooth
        isonormals(Xm, Ym, Zm, double(ctVol), h_metal);
        set(h_metal, 'EdgeColor','none', 'FaceColor', [0.80 0.20 0.20], 'FaceAlpha', 0.20);
        V0 = h_metal.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras.'; h_metal.Vertices = Vras(:,1:3);
    
        % Plot C and u (in RAS)
        plot3(C(1), C(2), C(3), 'ko', 'MarkerFaceColor','y', 'MarkerSize',8, ...
              'DisplayName','C (center)');
        quiver3(C(1), C(2), C(3), vecL*u(1), vecL*u(2), vecL*u(3), 0, ...
                'r', 'LineWidth', 2, 'DisplayName','u (axis)');
    
        xlabel('X (RAS)'); ylabel('Y (RAS)'); zlabel('Z (RAS)');
        camlight headlight; lighting gouraud; legend('Location','best');
    end


    % diagnostic plot
    % thrHU = 3000;         % threshold for “metal”
    % vizR  = 12;           % crop half-width in mm (e.g., 10–15 mm works well)
    % vecL  = 10;           % arrow length in mm
    % showSamples = true;   % scatter some core samples on the ring
    % visualize_pose_ct(C, u, rad_mm, len_mm, Avox2ras, ctVol, voxdims, thrHU, vizR, vecL, showSamples);

    % axial taper (Hann) to suppress fringe contamination
    w_ax = 0.5*(1 - cos(2*pi*( (t - t(1)) / (t(end)-t(1)+eps) )));  % same length as t
    w_ax    = 0.5 + 0.5*w_ax;   % now ranges 0.5..1

    Lc = len_mm/2;                     % half contact length
    Ls = insulating_spacer_length;     % full spacer length (per side) used in t_shell
    
    aL = -(Lc + Ls);   bL = -Lc;       % left cap interval
    aR =  +(Lc);       bR =  Lc + Ls;  % right cap interval
    
    w_shell = zeros(size(t_shell));    % zero in (-Lc, +Lc) by default
    
    % Left cap: 0.5 at aL/bL, 1 at midpoint
    base_w_value = 0.2;
    maskL = (t_shell >= aL) & (t_shell <= bL);
    sL = (t_shell(maskL) - aL) / max(bL - aL, eps);    % 0..1
    w_shell(maskL) = base_w_value + 0.5 * cos(pi*(sL - 0.5));   % 0.5 -> 1 -> 0.5
    
    % Right cap: symmetric
    maskR = (t_shell >= aR) & (t_shell <= bR);
    sR = (t_shell(maskR) - aR) / max(bR - aR, eps);    % 0..1
    w_shell(maskR) = base_w_value + 0.5 * cos(pi*(sR - 0.5));   % 0.5 -> 1 -> 0.5
    midMask = (t_shell > -Lc) & (t_shell < +Lc);
    w_shell(midMask) = base_w_value;

    % w_shell= 0.5*(1 - cos(2*pi*( (t_shell - t_shell(1)) /(t_shell(end)- t_shell(1) + eps) )));


    th = linspace(0, 2*pi, 24);   % a bit denser angularly, angular samping, 24 axial sample, with 10 numel sample

    % inside cylinder
    r_in = max(rad_mm, 0.5);
    n_r_in = ceil(max(r_in, step_mm)/step_mm);% radial samples inside (0..rad_mm)
    [TT_in, TH_in, RR_in] = cylinder_grid(t, th, 0.0, r_in, n_r_in);
    % represent the coordinates (in cylindrical space) of the core (inside) cylinder, where:
	% 	TT_in: axial positions (t) — like Z along the electrode
	% 	TH_in: angular positions (θ) — from 0 to 2π
	% 	RR_in: radial positions (r) — from 0 up to r_in
    % 
    Pin = cyl_points(C,u,v1,v2,TT_in,TH_in,RR_in);

    % the outside cylinder from 1.5 to 2.5 * rad_mm, for 
    Lc = len_mm/2;
    Ls = insulating_spacer_length;              % per your t_shell definition
    
    % Mid-annulus axial vector (within contact)
    t_mid = (-Lc):step_mm:(+Lc);                % subset of t_shell region
    % Caps axial vectors (outside contact, within spacer)
    t_capL = (-(Lc + Ls)):step_mm:(-Lc);
    t_capR = (+Lc):step_mm:(+(Lc + Ls));
    
    % ---- radial ranges ----
    % (1) Mid annulus: 1.5–2.5 × rad_mm (ensure outside the core)
    r_mid_lo = max(1.5*rad_mm, r_in + 0.2);   % keep a small safety margin
    r_mid_hi = 2.5 * rad_mm;
    
    % (2) Caps: full disk 0–2.5 × rad_mm
    r_cap_lo = 0.0;
    r_cap_hi = 2.5 * rad_mm;
    
    % Radial levels (voxel-aware)
    n_r_mid = max(3, ceil((r_mid_hi - r_mid_lo)/step_mm));
    n_r_cap = max(4, ceil((r_cap_hi - r_cap_lo)/step_mm));
    
    % ---- build shell samples (mid annulus) ----
    if ~isempty(t_mid)
        [TT_mid, TH_mid, RR_mid] = cylinder_grid(t_mid, th, r_mid_lo, r_mid_hi, n_r_mid);
    else
        TT_mid = []; TH_mid = []; RR_mid = [];
    end
    
    % ---- build shell samples (caps full disk) ----
    [TT_capL, TH_capL, RR_capL] = cylinder_grid(t_capL, th, r_cap_lo, r_cap_hi, n_r_cap);
    [TT_capR, TH_capR, RR_capR] = cylinder_grid(t_capR, th, r_cap_lo, r_cap_hi, n_r_cap);
    
    % ---- concatenate shell samples (mid + caps) ----
    TT_sh = [TT_mid; TT_capL; TT_capR];
    TH_sh = [TH_mid; TH_capL; TH_capR];
    RR_sh = [RR_mid; RR_capL; RR_capR];
    
    % 3D shell points
    Psh = cyl_points(C,u,v1,v2,TT_sh,TH_sh,RR_sh);

    % Sample CT (correct order inside helper)
    vin = sample_world_interp3(ctVol, Avox2ras, Pin);
    vsh = sample_world_interp3(ctVol, Avox2ras, Psh);

    % Apply axial taper weights
    % sampling along the axis t. The ends of the sampling range (near ±len/2 or the shell caps) are where artifacts are most likely:
	% partial‑volume & blur at the contact edges,
	% tiny mis‑registration of C,u shifts samples off the lead,
	% CT streaks/bloom broaden near metal transitions.
    % wv = repmat(w_ax(:), numel(th), 1);   % repeat along theta
    wv_in = repmat(w_ax(:), numel(th)*n_r_in, 1);
    idx_shell = interp1(t_shell, 1:numel(t_shell), TT_sh, 'nearest', 'extrap');
    idx_shell = max(1, min(numel(w_shell), idx_shell));
    wv_sh = w_shell(idx_shell).';
    
    vin_w = vin .* wv_in;
    vsh_w = vsh .* wv_sh;

    % pass weights scheck
    do_viz_weights = false;
    if do_viz_weights
        figure('Color','w'); hold on; grid on; box on;
        plot(t_shell, w_shell, 'LineWidth', 2);
    
        yl = ylim;
        % Shade left/right caps and middle region for clarity
        patch([aL bL bL aL],[yl(1) yl(1) yl(2) yl(2)], [0.85 0.95 1.0], ...
              'EdgeColor','none', 'FaceAlpha', 0.25);
        patch([aR bR bR aR],[yl(1) yl(1) yl(2) yl(2)], [0.85 0.95 1.0], ...
              'EdgeColor','none', 'FaceAlpha', 0.25);
        patch([-Lc +Lc +Lc -Lc],[yl(1) yl(1) yl(2) yl(2)], [0.98 0.88 0.88], ...
              'EdgeColor','none','FaceAlpha',0.25);
    
        % Reference lines
        yline(1.0,'k:'); yline(0.5,'k:'); xline(-Lc,'k--'); xline(+Lc,'k--');
        xline(aL,'b--'); xline(bL,'b--'); xline(aR,'b--'); xline(bR,'b--');
    
        xlabel('t (mm) along axis');
        ylabel('w\_shell');
        title('Shell axial weights with cap/middle regions');
        legend({'w\_shell','left cap','right cap','middle annulus'}, 'Location','best');
    end

    do_viz_sampling = false;
    if do_viz_sampling
        figure('Color','w'); hold on; axis equal vis3d; grid on; box on;
    
        % Optional: local metal iso‑crop for context (fast & light)
        try
            thrHU_viz = 2500; vizR = 12;
            Ainv = inv(Avox2ras); ijk0 = [C,1]*Ainv.'; i0=ijk0(2)+1; j0=ijk0(1)+1; k0=ijk0(3)+1;
            [m,n,p] = size(ctVol);
            hi = max(2, round(vizR / max(voxdims(1),eps)));
            hj = max(2, round(vizR / max(voxdims(2),eps)));
            hk = max(2, round(vizR / max(voxdims(3),eps)));
            Ir = max(1,floor(i0-hi)):min(m,ceil(i0+hi));
            Jr = max(1,floor(j0-hj)):min(n,ceil(j0+hj));
            Kr = max(1,floor(k0-hk)):min(p,ceil(k0+hk));
            BW = ctVol(Ir,Jr,Kr) > thrHU_viz;
            [Xc,Yc,Zc] = meshgrid(Jr,Ir,Kr);
            if any(BW(:))
                hIso = patch(isosurface(Xc,Yc,Zc,BW,0.5));
                isonormals(Xc,Yc,Zc,double(BW),hIso);
                set(hIso,'EdgeColor','none','FaceAlpha',0.18,'FaceColor',[0.7 0.7 0.9]);
                V0 = hIso.Vertices - 1;
                Vras = [V0, ones(size(V0,1),1)] * Avox2ras.';
                hIso.Vertices = Vras(:,1:3);
            end
        catch
            % safe to ignore viz errors
        end
    
        % Axis arrow
        q = 8; u_unit = u./max(norm(u),eps);
        quiver3(C(1),C(2),C(3), q*u_unit(1), q*u_unit(2), q*u_unit(3), 0, 'r', 'LineWidth',1.5);
    
        % Decimate for clarity
        dec_in = max(1, round(size(Pin,1)/2000));
        dec_sh = max(1, round(size(Psh,1)/3000));
    
        % Core points
        plot3(Pin(1:dec_in:end,1), Pin(1:dec_in:end,2), Pin(1:dec_in:end,3), ...
        '.', 'Color',[0.1 0.7 0.1], 'MarkerSize',8, 'DisplayName','Core');

        % Shell points colored by weight (wv_sh)
        idx = 1:dec_sh:numel(wv_sh);
        scatter3(Psh(idx,1), Psh(idx,2), Psh(idx,3), ...
                 14, wv_sh(idx), 'filled', 'Marker','^', 'MarkerEdgeColor','k', ...
                 'DisplayName','Shell (weighted)');
        
        % Colormap / colorbar for weights
        colormap(parula);
        cb = colorbar; ylabel(cb,'Shell axial weight');
        clim([min(wv_sh) max(wv_sh)]);  % or [0 1] if you prefer fixed scaling
        
        xlabel('X (RAS)'); ylabel('Y (RAS)'); zlabel('Z (RAS)');
        title('Shell samples colored by axial weight (core in green)');
        camlight headlight; lighting gouraud;
        legend('Location','best');
    end


    % check where the core is and where the shell is
    
    % Robust contrast: high inside vs median shell, normalized by robust scale
    in_p50 = prctile(vin_w, 50);
    sh_p50 = prctile(vsh_w, 50);
    sh_iqr = iqr(vsh_w(~isnan(vsh_w)));
    sh_sd  = max(sh_iqr/1.349, 1);        % avoid divide-by-small
    mu     = (in_p50 - sh_p50) / sh_sd;
end



function vals = sample_world_interp3(vol, Avox2ras, P)
% no need switch u  v w 
    Ainv = inv(Avox2ras);
    uvw1 = [P, ones(size(P,1),1)] * Ainv.';  % -> voxel (0-based)
    u = uvw1(:,1)+1; v = uvw1(:,2)+1; w = uvw1(:,3)+1;
    vals = interp3(vol, u, v, w, 'linear', NaN);
end

function [TT,TH,RR] = cylinder_grid(t, th, r_lo, r_hi, n_r)
% Quasi-uniform sampling over cylindrical cross-section using sqrt spacing.
    [TT,TH] = ndgrid(t, th);
    r2_lo = max(0, r_lo^2); r2_hi = max(r2_lo + eps, r_hi^2);
    r_levels = sqrt(linspace(r2_lo, r2_hi, n_r));      % 1 x n_r
    TT = repmat(TT, [1, 1, n_r]);
    TH = repmat(TH, [1, 1, n_r]);
    RR = repmat(reshape(r_levels,1,1,[]), [numel(t), numel(th), 1]);
    TT = TT(:); TH = TH(:); RR = RR(:);
end

function P = cyl_points(C,u,v1,v2,TT,TH,RR)
    n = numel(TT);
    P = repmat(C, n,1) + TT.*repmat(u,n,1) + RR.*cos(TH).*repmat(v1,n,1) + RR.*sin(TH).*repmat(v2,n,1);
end



    % old script
    % n_r_sh = 4;
    % [TT_sh, TH_sh, RR_sh] = cylinder_grid(t_shell, th, r_in, r_out, n_r_sh);
    % 
    % 
    % [TT,TH] = ndgrid(t, th);
    % TTv = TT(:); THv = TH(:);            % Nx1
    % N = numel(TTv);
    % 
    % % Expand direction rows to Nx3
    % C1  = repmat(C,  N,1);
    % u1  = repmat(u,  N,1);
    % v1r = repmat(v1, N,1);
    % v2r = repmat(v2, N,1);
    % 
    % % ----- inside cylinder points -----
    % 
    % rin_cos = rin .* cos(THv);
    % rin_sin = rin .* sin(THv);
    % 
    % Uterm = TTv .* u1;
    % R1    = rin_cos .* v1r;
    % R2    = rin_sin .* v2r;
    % Pin   = C1 + Uterm + R1 + R2;        % Nx3
    % 
    % % ----- shell annulus (use midpoint between inner & outer radii) -----
    % r_in_shell  = 1.5 * rad_mm;            % move shell farther out
    % r_out_shell = 3.0*rad_mm;
    % r_mid       = 0.5*(r_in_shell + r_out_shell);
    % 
    % rmc = r_mid .* cos(THv);
    % rms = r_mid .* sin(THv);
    % R1s = rmc .* v1r;
    % R2s = rms .* v2r;
    % Psh = C1 + Uterm + R1s + R2s;
% showRingOnly = false;   % show only central axial ring to keep it clean
% 
%     % --- crop a local subvolume around C (convert mm->vox per axis) ---
%     if showRingOnly
%         Ainv = inv(Avox2ras);
%         ijk0 = [C,1]*Ainv.';      % 0-based [i j k]
%         i0 = ijk0(1)+1; j0 = ijk0(2)+1; k0 = ijk0(3)+1;
% 
%         [m,n,p] = size(ctVol);
%         hi = max(2, round(vizR / max(voxdims(1),eps)));
%         hj = max(2, round(vizR / max(voxdims(2),eps)));
%         hk = max(2, round(vizR / max(voxdims(3),eps)));
% 
%         Ir = max(1,floor(i0-hi)) : min(m,ceil(i0+hi));
%         Jr = max(1,floor(j0-hj)) : min(n,ceil(j0+hj));
%         Kr = max(1,floor(k0-hk)) : min(p,ceil(k0+hk));
% 
%         subMask = ctVol(Ir, Jr, Kr) > thrHU;
% 
%         % --- isosurface on crop, using ABSOLUTE 1-based indices ---
%         [Xc, Yc, Zc] = meshgrid(Jr, Ir, Kr);   % X=j, Y=i, Z=k  (ABSOLUTE, 1-based)
% 
%         figure('Color','w'); hold on; axis equal vis3d; grid on; box on;
%         title('Zoom: metal isosurface, C/u, core & shell samples');
% 
%         hIso = patch(isosurface(Xc, Yc, Zc, subMask, 0.5));
%         isonormals(Xc, Yc, Zc, double(subMask), hIso);
%         set(hIso,'EdgeColor','none','FaceAlpha',0.20,'FaceColor',[0.7 0.7 0.9]);
% 
%         % Vertices are [X Y Z] = [j i k] in 1-based absolute voxel coords
%         V1 = hIso.Vertices;
%         V0 = V1 - 1;  % -> absolute 0-based [j0 i0 k0], consistent with full-volume path
% 
%         % Map voxel(0-based) -> RAS (no axis swap)
%         Vras = [V0, ones(size(V0,1),1)] * Avox2ras.';
%         hIso.Vertices = Vras(:,1:3);
% 
%         % [m,n,p] = size(ctVol);
%         % [Xm, Ym, Zm] = meshgrid(1:n, 1:m, 1:p);
%         % BW = ctVol > thrHU;
%         % h_metal = patch(isosurface(Xm, Ym, Zm, BW, 0.5));
%         % % make 3d surface look smooth
%         % isonormals(Xm, Ym, Zm, double(ctVol), h_metal);
%         % set(h_metal, 'EdgeColor','none', 'FaceColor', [0.80 0.20 0.20], 'FaceAlpha', 0.20);
%         % V0 = h_metal.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras.'; h_metal.Vertices = Vras(:,1:3);
% 
%         % Plot C and u
%         plot3(C(1),C(2),C(3),'ko','MarkerFaceColor','y','MarkerSize',8);
%         quiver3(C(1),C(2),C(3),vecL*u(1),vecL*u(2),vecL*u(3),0,'r','LineWidth',2);
%     end
% 
%     % --- choose which sample points to draw to reduce clutter ---
%     if showRingOnly
%         % central axial ring (t ≈ 0)
%         midMask = abs(TTv) <= (max(voxdims)/2);   % ± half a voxel along axis
%         Pin_show = Pin(midMask,:);
%         Psh_show = Psh(midMask,:);
%         vin_show = vin(midMask);
%         vsh_show = vsh(midMask);
%     else
%         % decimate all samples
%         dec = max(1, round(numel(vin)/400));   % ~400 points max
%         idx = 1:dec:numel(vin);
%         Pin_show = Pin(idx,:);  Psh_show = Psh(idx,:);
%         vin_show = vin(idx);    vsh_show = vsh(idx);
%     end
% 
%     if showRingOnly
%     % Common color scaling from pooled samples
%         pool = [vin_show(:); vsh_show(:)];
%         cmin = prctile(pool, 5);  cmax = prctile(pool, 95);
%         if ~isfinite(cmin) || ~isfinite(cmax) || cmin>=cmax, cmin = min(pool); cmax=max(pool); end
% 
%         % Core samples (magenta-ish markers)
%         sc1 = scatter3(Pin_show(:,1), Pin_show(:,2), Pin_show(:,3), ...
%                        18, vin_show, 'filled', 'MarkerEdgeColor','k', 'DisplayName','core samples');
%         % Shell samples (cyan-ish markers)
%         sc2 = scatter3(Psh_show(:,1), Psh_show(:,2), Psh_show(:,3), ...
%                        18, vsh_show, 'filled', 'MarkerEdgeColor','k', 'Marker','^', 'DisplayName','shell samples');
% 
%         colormap(parula); caxis([cmin cmax]); cb = colorbar; ylabel(cb,'HU');
% 
% 
% 
%         xlabel('X (RAS)'); ylabel('Y (RAS)'); zlabel('Z (RAS)');
%         camlight headlight; lighting gouraud;
%         legend('Location','best');
%         view(3);
%     end



