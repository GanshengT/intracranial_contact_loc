function mu = score_pose_contrast(C,u,rad_mm,len_mm,Avox2ras,ctVol,voxdims)
% cylinder coordinates: (t,r,\theta)=C+t\,u+r(\cos\theta\,v_1+\sin\theta\,v_2).
    % ---- input hygiene ----
    C = reshape(double(C),1,3);
    u = reshape(double(u),1,3);
    u = u./max(norm(u),eps);

    [v1,v2] = orthobasis_row(u);   % 1x3 each

    % Sampling density tied to native voxel sizes
    step_mm = max(min(voxdims)/2, 0.2);
    t  = -len_mm/3*2 : step_mm : +len_mm/3*2; % introduce slightly overlap

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

    th = linspace(0, 2*pi, 24);   % a bit denser angularly, angular samping, 24 axial sample, with 10 numel sample

    [TT,TH] = ndgrid(t, th);
    TTv = TT(:); THv = TH(:);            % Nx1
    N = numel(TTv);

    % Expand direction rows to Nx3
    C1  = repmat(C,  N,1);
    u1  = repmat(u,  N,1);
    v1r = repmat(v1, N,1);
    v2r = repmat(v2, N,1);

    % ----- inside cylinder points -----
    rin = 1.5 * rad_mm;
    rin_cos = rin .* cos(THv);
    rin_sin = rin .* sin(THv);

    Uterm = TTv .* u1;
    R1    = rin_cos .* v1r;
    R2    = rin_sin .* v2r;
    Pin   = C1 + Uterm + R1 + R2;        % Nx3

    % ----- shell annulus (use midpoint between inner & outer radii) -----
    r_in_shell  = 2 * rad_mm;            % move shell farther out
    r_out_shell = 3.0*rad_mm;
    r_mid       = 0.5*(r_in_shell + r_out_shell);

    rmc = r_mid .* cos(THv);
    rms = r_mid .* sin(THv);
    R1s = rmc .* v1r;
    R2s = rms .* v2r;
    Psh = C1 + Uterm + R1s + R2s;

    % Sample CT (correct order inside helper)
    vin = sample_world_interp3(ctVol, Avox2ras, Pin);
    vsh = sample_world_interp3(ctVol, Avox2ras, Psh);

    % Apply axial taper weights
    wv = repmat(w_ax(:), numel(th), 1);   % repeat along theta
    vin_w = vin .* wv;
    vsh_w = vsh .* wv;

    % check where the core is and where the shell is
    thrHU  = 3000;   % metal threshold (HU)
    vizR   = 60;     % half-width crop (mm) around C for the isosurface
    vecL   = 10;     % arrow length for axis
    showRingOnly = false;   % show only central axial ring to keep it clean
    
    % --- crop a local subvolume around C (convert mm->vox per axis) ---
    if showRingOnly
        Ainv = inv(Avox2ras);
        ijk0 = [C,1]*Ainv.';      % 0-based [i j k]
        i0 = ijk0(1)+1; j0 = ijk0(2)+1; k0 = ijk0(3)+1;
        
        [m,n,p] = size(ctVol);
        hi = max(2, round(vizR / max(voxdims(1),eps)));
        hj = max(2, round(vizR / max(voxdims(2),eps)));
        hk = max(2, round(vizR / max(voxdims(3),eps)));
        
        Ir = max(1,floor(i0-hi)) : min(m,ceil(i0+hi));
        Jr = max(1,floor(j0-hj)) : min(n,ceil(j0+hj));
        Kr = max(1,floor(k0-hk)) : min(p,ceil(k0+hk));
        
        subMask = ctVol(Ir, Jr, Kr) > thrHU;
        
        % --- isosurface on crop, using ABSOLUTE 1-based indices ---
        [Xc, Yc, Zc] = meshgrid(Jr, Ir, Kr);   % X=j, Y=i, Z=k  (ABSOLUTE, 1-based)
        
        figure('Color','w'); hold on; axis equal vis3d; grid on; box on;
        title('Zoom: metal isosurface, C/u, core & shell samples');
        
        hIso = patch(isosurface(Xc, Yc, Zc, subMask, 0.5));
        isonormals(Xc, Yc, Zc, double(subMask), hIso);
        set(hIso,'EdgeColor','none','FaceAlpha',0.20,'FaceColor',[0.7 0.7 0.9]);
        
        % Vertices are [X Y Z] = [j i k] in 1-based absolute voxel coords
        V1 = hIso.Vertices;
        V0 = V1 - 1;  % -> absolute 0-based [j0 i0 k0], consistent with full-volume path
        
        % Map voxel(0-based) -> RAS (no axis swap)
        Vras = [V0, ones(size(V0,1),1)] * Avox2ras.';
        hIso.Vertices = Vras(:,1:3);
    
        % [m,n,p] = size(ctVol);
        % [Xm, Ym, Zm] = meshgrid(1:n, 1:m, 1:p);
        % BW = ctVol > thrHU;
        % h_metal = patch(isosurface(Xm, Ym, Zm, BW, 0.5));
        % % make 3d surface look smooth
        % isonormals(Xm, Ym, Zm, double(ctVol), h_metal);
        % set(h_metal, 'EdgeColor','none', 'FaceColor', [0.80 0.20 0.20], 'FaceAlpha', 0.20);
        % V0 = h_metal.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras.'; h_metal.Vertices = Vras(:,1:3);
        
        % Plot C and u
        plot3(C(1),C(2),C(3),'ko','MarkerFaceColor','y','MarkerSize',8);
        quiver3(C(1),C(2),C(3),vecL*u(1),vecL*u(2),vecL*u(3),0,'r','LineWidth',2);
    end
        
    % --- choose which sample points to draw to reduce clutter ---
    if showRingOnly
        % central axial ring (t ≈ 0)
        midMask = abs(TTv) <= (max(voxdims)/2);   % ± half a voxel along axis
        Pin_show = Pin(midMask,:);
        Psh_show = Psh(midMask,:);
        vin_show = vin(midMask);
        vsh_show = vsh(midMask);
    else
        % decimate all samples
        dec = max(1, round(numel(vin)/400));   % ~400 points max
        idx = 1:dec:numel(vin);
        Pin_show = Pin(idx,:);  Psh_show = Psh(idx,:);
        vin_show = vin(idx);    vsh_show = vsh(idx);
    end
    if showRingOnly
    % Common color scaling from pooled samples
        pool = [vin_show(:); vsh_show(:)];
        cmin = prctile(pool, 5);  cmax = prctile(pool, 95);
        if ~isfinite(cmin) || ~isfinite(cmax) || cmin>=cmax, cmin = min(pool); cmax=max(pool); end
        
        % Core samples (magenta-ish markers)
        sc1 = scatter3(Pin_show(:,1), Pin_show(:,2), Pin_show(:,3), ...
                       18, vin_show, 'filled', 'MarkerEdgeColor','k', 'DisplayName','core samples');
        % Shell samples (cyan-ish markers)
        sc2 = scatter3(Psh_show(:,1), Psh_show(:,2), Psh_show(:,3), ...
                       18, vsh_show, 'filled', 'MarkerEdgeColor','k', 'Marker','^', 'DisplayName','shell samples');
        
        colormap(parula); caxis([cmin cmax]); cb = colorbar; ylabel(cb,'HU');
        
    
        
        xlabel('X (RAS)'); ylabel('Y (RAS)'); zlabel('Z (RAS)');
        camlight headlight; lighting gouraud;
        legend('Location','best');
        view(3);
    end


    % Robust contrast: high inside vs median shell, normalized by robust scale
    in_p80 = prctile(vin_w, 50);
    sh_p50 = prctile(vsh_w, 50);
    sh_iqr = iqr(vsh_w(~isnan(vsh_w)));
    sh_sd  = max(sh_iqr/1.349, 1);        % avoid divide-by-small
    mu     = (in_p80 - sh_p50) / sh_sd;
end



function vals = sample_world_interp3(vol, Avox2ras, P)
% no need switch u  v w 
    Ainv = inv(Avox2ras);
    uvw1 = [P, ones(size(P,1),1)] * Ainv.';  % -> voxel (0-based)
    u = uvw1(:,1)+1; v = uvw1(:,2)+1; w = uvw1(:,3)+1;
    vals = interp3(vol, u, v, w, 'linear', NaN);
end