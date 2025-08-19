function [G, info] = grow_from_seed( ...
    seed_ijk, ctVol, ROI, metal_roi_lo, Avox2ras0, ...
    start_traj, U_in, max_voxels, vx, varargin)
% GROW_FROM_SEED
% g.tan@wustl.edu
% This function will grow a connected bolb from a seed
% one input is the direction of the planned line
% from the current growing bolb, we will search for a cloest voxel that is
% 1 in metal_roi_lo, then we will decide whether we will include it, add to
% visited set. 
% 
% Acceptance:
%   - If voxel in metal_roi_lo => accept (anchor to metal).
%   - Else accept if score > score_min (defaults to >0).
%
% Priority score for candidate v:
%   score(v) = w_fill * fill%  +  w_axis * (1 - r(v)/r_outer)  +  w_fwd * fwd_term
% where fill% of a fitted t-bin if we include v
% based on the current growing bolb, we will get the centroid. we also have
% the planned trajector direction (we will allow direction around the
% planned traj), then we will get a tube along the planned trajectory
% direction (we will later search direction, such that fill% is maximized),
% the tube will be fitted based on centroid, height and the axis direction

% question here, do we need w_axis and w_fwd? Alternatively, we check if
% average distance between the candidate voxel and existing bolb is larger
% than 3 SD of the inter-voxel distance within the existing bolb. this can
% be added to the composite score

% KNearest is set to 1, 1 at a time

 % ---------- options ----------
    ip = inputParser; ip.CaseSensitive=false;
    addParameter(ip,'UseROI',false,@(x)islogical(x)||ismember(x,[0 1]));
    addParameter(ip,'AngleDeg',4,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=45);
    addParameter(ip,'DirGrid',[5 8],@(x)isnumeric(x)&&numel(x)==2);
    addParameter(ip,'TPadMM',5.0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'Weights',[1 0.1],@(x)isnumeric(x)&&numel(x)==2);
    addParameter(ip,'SigmaK',3.0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'ScoreMin',-0.2,@(x)isnumeric(x)&&isscalar(x));
    addParameter(ip,'FillMin',0.1,@(x)isnumeric(x)&&isscalar(x));
    addParameter(ip,'fill_drop_tol',0.2,@(x)isnumeric(x)&&isscalar(x));
    addParameter(ip,'KNearest',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'Trace',true,@(x)islogical(x)||ismember(x,[0 1]));
    addParameter(ip,'TraceN',300,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'KInit',8,@(x)isnumeric(x)&&isscalar(x)&&x>=0);     % how many initial nearest voxels
    addParameter(ip,'InitUseROI',false,@(x)islogical(x)||ismember(x,[0 1])); % constrain KInit to ROI?
    parse(ip,varargin{:});
    UseROI   = logical(ip.Results.UseROI);
    AngleDeg = ip.Results.AngleDeg;
    DirGrid  = ip.Results.DirGrid;
    TPadMM   = ip.Results.TPadMM;
    W        = ip.Results.Weights(:).';
    ScoreMin = ip.Results.ScoreMin;
    KNearest = ip.Results.KNearest;
    FillMin   = ip.Results.FillMin;
    fill_drop_tol   = ip.Results.fill_drop_tol;
    KInit     = ip.Results.KInit;
    InitUseROI= logical(ip.Results.InitUseROI);
    Trace= logical(ip.Results.Trace);


    log_cap = 256;
    logs = repmat(struct( ...
        'step',0, 'ijk',[NaN NaN NaN], 'ras',[NaN NaN NaN], ...
        'decision',"",'reason',"", 'dFill',NaN, 'fill_new',NaN, 'cur_fill',NaN, ...
        'score',NaN, 'dist_pen',NaN, 'dmin',NaN, 'dave',NaN, ...
        't0',NaN, 't1',NaN, 'r_mm',NaN, 'ctr',[NaN NaN NaN], 'dstar',[NaN NaN NaN], ...
        'time',NaN, 'note',""), 1, log_cap);
    log_n = 0;
    

    % ---------- checks ----------
    ROI          = logical(ROI);
    metal_roi_lo = logical(metal_roi_lo);
    assert(isequal(size(ctVol), size(ROI), size(metal_roi_lo)), 'size mismatch');
    [M,N,P] = size(ctVol);

    I = seed_ijk(1); J = seed_ijk(2); K = seed_ijk(3);
    assert(I>=1 && I<=M && J>=1 && J<=N && K>=1 && K<=P, 'seed_ijk OOB');

    U = U_in(:).'; U = U / max(norm(U),eps);

    % ---------- precompute voxel centers & t projection ----------
    [Xras, t_field] = precompute_ras_t(Avox2ras0, start_traj, U, M,N,P);

    % ---------- init ----------
    G = false(M,N,P);
    visited = false(M,N,P);
    rejected = false(M,N,P);

    % seed
    G(I,J,K) = true; visited(I,J,K) = true;

    % candidate base mask (metal AND optional ROI)
    cand_base = metal_roi_lo;
    if UseROI, cand_base = cand_base & ROI; end

    % diagnostics
    rc = struct('reject_low_score',0,'empty_candidates',0,'cap_reached',0);
    % samp = struct('low_score',zeros(0,3));
    % KMAX = 200;


    grown = 1;
    stop_cause = 'pq_empty';

    % --- per-slab fill diagnostics (based on metal availability) ---
    % [edges, cap_per_bin] = summarize_fill_metal(cand_base, t_field, TBinMM);
    % t_field is a scalar field with the same shape as the CT, giving each
    % voxel’s longitudinal coordinate (in mm) along the planned trajectory
    % axis.
    % acc_per_bin = zeros(size(cap_per_bin));
    if KInit > 0
        % candidate pool for prefill: metal (and optional ROI), not visited/G
        cand_init = metal_roi_lo & ~visited & ~G;
        if InitUseROI, cand_init = cand_init & ROI; end
    
        % KInit nearest to the seed point (RAS distance)
        seed_ras = squeeze(Xras(I,J,K,:)).';
        [iiC,jjC,kkC, dC] = k_nearest_to_point(cand_init, Xras, seed_ras, KInit);
    
        % accept them (unconditionally) and update diagnostics
        for q = 1:numel(iiC)
            ii = iiC(q); jj = jjC(q); kk = kkC(q);
            G(ii,jj,kk) = true; visited(ii,jj,kk) = true;
        end
        grown = grown + numel(iiC);
        if grown >= max_voxels
            rc.cap_reached = rc.cap_reached + 1; stop_cause = 'cap_reached';
            goto_main_loop = false;
        else
            goto_main_loop = true;
        end
    else
        goto_main_loop = true;
    end
    
    if ~goto_main_loop
        % bail early if cap hit during prefill
        % (info pack & viz happen later as already implemented)
        return;
    end

    % diagnose initial bolb
    [m,n,p] = size(ctVol); [Xm,Ym,Zm] = meshgrid(1:n,1:m,1:p);
    % figure('Color','w'); hold on;
    % % grown surface
    % if any(G(:))
    %     hg = patch(isosurface(Xm,Ym,Zm,G,0.5));
    %     set(hg,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.5);
    %     V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
    %     hg.Vertices = Vras(:,1:3);
    % end
    % hg = patch(isosurface(Xm,Ym,Zm,metal_roi_lo,0.5));
    % set(hg,'EdgeColor','none','FaceColor',[1 0.7 0],'FaceAlpha',0.3);
    % V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
    % hg.Vertices = Vras(:,1:3);

    % ---------- main loop ----------
    tube_change = true;
    nn_change   = true;

    ctr = [NaN NaN NaN]; dstar = [NaN NaN NaN];
    t0 = NaN; t1 = NaN; r_mm = NaN; cur_fill = NaN;
    d_vec_all = []; muD = NaN; sigD = 1;   % guards
    while true
        % (1) Tube fitting around current blob: axis direction d*, center at blob centroid
        % ctr = blob_centroid_ras(G, Xras); % your existing helper
        % 
        % [dstar, cur_fill, t0, t1, r_mm] = fit_axis_by_fill( ...
        %     U, ctr, TPadMM, G, Xras, vx, ...
        %     'AngleDeg', AngleDeg, 'RingK', 12, 'IncludeCenter', true, ...
        %     'RGuardMM', 0.2, 'TGuardMM', 0.5, 'RobustPct', 98);
        % 
        % [t0, t1, r_mm, stats] = enclosing_tube_params(G, Xras, ctr, dstar);

        if tube_change
            ctr = blob_centroid_ras(G, Xras);
            [dstar, ~, t0, t1, r_mm] = fit_axis_by_fill( ...
                U, ctr, TPadMM, G, Xras, vx, ...
                'AngleDeg', AngleDeg, 'RingK', 12, 'IncludeCenter', true, ...
                'RGuardMM', 0.2, 'TGuardMM', 0.5, 'RobustPct', 98);
    
            % tighten with enclosing tube
            [t0, t1, r_mm, stats] = enclosing_tube_params(G, Xras, ctr, dstar);
    
            % fill% using current tube
            V_tube_mm3 = max(1e-9, pi * (r_mm^2) * max(0, (t1 - t0)));
            cur_fill = safe_ratio(nnz(G) * prod(vx), V_tube_mm3);
    
            tube_change = false;   % now up-to-date
        end


        % build its mask (fast if you pass Domain = cand_base | G)
        % Domain = (cand_base | G);
        % tube_mask = tube_mask_from_params(ctr, dstar, t0, t1, r_enclose_mm, Xras, Domain);

        % % diagnose: 
        % figure('Color','w'); hold on;
        % % Blob surface (G) in RAS
        % if any(G(:))
        %     [m,n,p] = size(G);
        %     [Xm,Ym,Zm] = meshgrid(1:n,1:m,1:p);
        %     hg = patch(isosurface(Xm,Ym,Zm,G,0.5));
        %     set(hg,'EdgeColor','none','FaceColor',[1 0.7 0],'FaceAlpha',0.35);
        %     V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
        %     hg.Vertices = Vras(:,1:3);
        % end
        % 
        % % Planned line (span length based on voxel bbox extent)
        % U0 = U_in(:).'/max(norm(U_in),eps);
        % % quick length: use G extent if available; else a default span
        % if any(G(:))
        %     [ii,jj,kk] = ind2sub(size(G), find(G));
        %     PG = [ squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*1))), ...
        %            squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*2))), ...
        %            squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*3))) ];
        %     tG = (PG - start_traj) * U0(:);
        %     L0 = [min(tG) max(tG)] + [-TPadMM TPadMM];
        % else
        %     L0 = [-10 10]; % mm fallback
        % end
        % P0 = start_traj + L0(1)*U0; P1 = start_traj + L0(2)*U0;
        % plot3([P0(1) P1(1)],[P0(2) P1(2)],[P0(3) P1(3)],'-','LineWidth',2,'Color',[0 0.45 0.74]);
        % if any(G(:))
        %     ctrPlot = mean(PG,1);
        %     scatter3(ctrPlot(1),ctrPlot(2),ctrPlot(3),80,[0.85 0.33 0.10],'filled','MarkerEdgeColor','k');
        % end
        % 
        % % hg = patch(isosurface(Xm,Ym,Zm,tube_mask,0.5));
        % % set(hg,'EdgeColor','none','FaceColor',[0 1 0],'FaceAlpha',0.35);
        % % V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
        % % hg.Vertices = Vras(:,1:3);
        % EPS_NUM = 1e-9;
        % d = dstar(:).'/max(norm(dstar),EPS_NUM);     % unit axis
        % 
        % % build an orthonormal frame {d, n1, n2}
        % tmp = [1 0 0]; if abs(dot(tmp,d))>0.9, tmp = [0 1 0]; end
        % n1 = tmp - dot(tmp,d)*d; n1 = n1 / max(norm(n1),EPS_NUM);
        % n2 = cross(d, n1);       n2 = n2 / max(norm(n2),EPS_NUM);
        % 
        % % longitudinal sampling (denser if longer tube)
        % L  = max(t1 - t0, 1e-3);
        % nTau = max(24, ceil(L/3));     % ~3 mm spacing along axis
        % nTh  = 48;                      % circular resolution
        % 
        % tau   = linspace(t0, t1, nTau);
        % theta = linspace(0, 2*pi, nTh);
        % [TH, TAU] = meshgrid(theta, tau);
        % 
        % % parametric tube surface in RAS
        % Xt = ctr(1) + TAU.*d(1) + r_mm*(cos(TH)*n1(1) + sin(TH)*n2(1));
        % Yt = ctr(2) + TAU.*d(2) + r_mm*(cos(TH)*n1(2) + sin(TH)*n2(2));
        % Zt = ctr(3) + TAU.*d(3) + r_mm*(cos(TH)*n1(3) + sin(TH)*n2(3));
        % 
        % hs_tube = surf(Xt, Yt, Zt, ...
        %     'FaceAlpha', 0.10, 'EdgeColor', 'none', 'FaceColor', [0 0.4 1]);
        % 
        % % (optional) end-caps to make the tube extent obvious
        % capTh = linspace(0, 2*pi, 64);
        % C0 = ctr + t0*d + r_mm*(cos(capTh)'*n1 + sin(capTh)'*n2);
        % C1 = ctr + t1*d + r_mm*(cos(capTh)'*n1 + sin(capTh)'*n2);
        % patch('XData',C0(:,1),'YData',C0(:,2),'ZData',C0(:,3), ...
        %       'FaceColor',[0 0.4 1],'FaceAlpha',0.10,'EdgeColor','none');
        % patch('XData',C1(:,1),'YData',C1(:,2),'ZData',C1(:,3), ...
        %       'FaceColor',[0 0.4 1],'FaceAlpha',0.10,'EdgeColor','none');
        % 
        % axis equal vis3d; grid on; box on; camlight headlight; lighting gouraud;
        % xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
        % title('Planned trajectory (blue), fitted tube (blue translucent), grown blob (orange)');
        % view(135,20);



        % current fill% for reporting
        % V_tube_mm3  = max(1e-9, pi * (r_mm^2) * max(0, (t1 - t0)));  
        % cur_fill = safe_ratio(nnz(G)* prod(vx), V_tube_mm3);

        % (2) choose K nearest **metal** candidates to the blob (not visited/rejected/in-G)
        eligible = cand_base & ~visited & ~rejected & ~G;
        if ~any(eligible(:))
            stop_cause = 'no candidates'; break;
        end
        [ic,jc,kc, dmin_vec, davg_vec] = k_nearest_to_blob( ...
             G, eligible, Xras, KNearest, 'SampleN', 2000, 'AvgMode', 'trimmed', 'TrimPct', 10);

        cand_list = [ic(:), jc(:), kc(:)];
        % nn stat change, like when we append a new candidate, recompute
        if nn_change 
            [muD, sigD, stats, d_vec_all] = nn_stats_blob(G, Xras);
            if sigD <= 1e-9, sigD = 1; end              
            nn_change = false;% guard
        end
        % Dthr = max(1e-6, muD + Ksigma*sigD);                    % soft cap for normalization
        
        % Thresholds and weights (tune)
        
        dist_max_mm     = 3.0;         % hard distance cap (mm) for a single jump
        alpha_fill      = W(1);         % weight for Δfill
        beta_dist       = W(2);         % weight for distance penalty
        
        % bestScore = -inf; bestIdx = []; bestParts = [0 0 0];  % [dFill, dMin, z]
        % for kk2 = 1:size(cand_list,1)
        kk2 = 1;
        ii = cand_list(kk2,1); jj = cand_list(kk2,2); kk = cand_list(kk2,3);
        dmin = dmin_vec(kk2);
        dave = davg_vec(kk2);        % robust soft metric

        % Hard guards first (optional but cheap)
        if dmin > dist_max_mm
            rejected(ii,jj,kk) = true;
            continue;                                       % too large a jump
        end
    
        [fill_new, dstar_tmp, ctr_tmp] = eval_fill_after_add( ...
            G, [ii jj kk], vx, U, TPadMM, ...
            cand_base, Xras, Avox2ras0, ...
            'AngleDeg', AngleDeg, 'RingK', 12, 'IncludeCenter', true, ...
            'RGuardMM', 0.2, 'TGuardMM', 0.5, 'RobustPct', 98);
        
        % V_tube_mm3  = max(1e-9, pi * (r_mm^2) * max(0, (t1 - t0)));  
        % fill_new = safe_ratio(nnz(G)* prod(vx), V_tube_mm3);

        dFill = fill_new - cur_fill;  
    
        % Soft distance penalty via z-score relative to current blob spacing
        if sigD <= 1e-9, sigD = 1; end

        if numel(d_vec_all) ==0 || ~isfinite(dave)
            dist_pen = 0;                               % safe fallback
        else
            % Empirical CDF: fraction of blob distances ≤ dave, in [0,1]
            dist_pen = sum(d_vec_all <= dave) / numel(d_vec_all);  
            % (If you prefer strict percentile, use < instead of ≤.)
        end
        score = alpha_fill * dFill - beta_dist * dist_pen;
        
        % log decision
        cand_ras = squeeze(Xras(ii,jj,kk,:)).';
        decision = "pending";
        % (B) OPTIONAL pre-compute "new tube" for diagnostics only
        t0_new = NaN; t1_new = NaN; r_new = NaN; ctr_new = ctr; dstar_new = dstar;
        if Trace
            Gtmp = G; Gtmp(ii,jj,kk) = true;
            try
                [t0_new, t1_new, r_new] = enclosing_tube_params(Gtmp, Xras, ctr_tmp, dstar_tmp);
                ctr_new   = ctr_tmp;
                dstar_new = dstar_tmp;
            catch
                % fall back silently if enclosing_tube_params fails
            end
        end

        % diagnostic plot
        
        % tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        % 
        % % (1) Current blob + CURRENT tube
        % nexttile; hold on;
        % [m_,n_,p_] = size(G); [Xm_,Ym_,Zm_] = meshgrid(1:n_,1:m_,1:p_);
        % if any(G(:))
        %     hg = patch(isosurface(Xm_,Ym_,Zm_,G,0.5));
        %     set(hg,'EdgeColor','none','FaceColor',[1 0.6 0.1],'FaceAlpha',0.35);
        %     V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
        %     hg.Vertices = Vras(:,1:3);
        % end
        % % draw current tube
        % local_draw_tube(ctr, dstar, t0, t1, r_mm);
        % % planned line
        % pad_mm = 1.0;                          % adjust as you like
        % U0 = U_in(:).'/max(norm(U_in),eps);    % unit planned direction
        % 
        % if any(G(:))
        %     tvals = t_field(G);                % mm along planned axis for blob voxels
        %     t0p = min(tvals) - pad_mm;         % padded min
        %     t1p = max(tvals) + pad_mm;         % padded max
        % else
        %     % fallback if blob empty: use current tube span (also padded)
        %     t0p = (t0 - pad_mm);
        %     t1p = (t1 + pad_mm);
        % end
        % 
        % P0 = start_traj + t0p*U0;
        % P1 = start_traj + t1p*U0;
        % 
        % plot3([P0(1) P1(1)], [P0(2) P1(2)], [P0(3) P1(3)], 'k--', 'LineWidth', 1.0);
        % scatter3(cand_ras(1), cand_ras(2), cand_ras(3), 36, 'r', 'filled');
        % axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
        % xlabel X; ylabel Y; zlabel Z;
        % title(sprintf('Step %d — Current blob + tube', grown));
        % 
        % % (2) Inter-voxel distance distribution + candidate percentile
        % nexttile; 
        % if ~isempty(d_vec_all)
        %     histogram(d_vec_all, 'BinMethod','fd'); hold on;
        %     yl = ylim;
        %     plot([dave dave], yl,'-','LineWidth',1.5);
        %     text(dave, yl(2)*0.95, sprintf('  d_{ave}=%.2f mm', dave), 'VerticalAlignment','top');
        %     title('Blob inter-voxel distances'); xlabel('mm'); ylabel('count');
        % else
        %     text(0.5,0.5,'No distance stats (tiny blob)','HorizontalAlignment','center');
        %     axis off;
        % end
        % 
        % % (3) Hypothetical NEW tube (blob + candidate)
        % nexttile; hold on;
        % if any(G(:))
        %     % Create a 3D logical array for the single voxel
        %     single_voxel = false(M, N, P);  % Assuming M, N, P are the dimensions
        %     single_voxel(ii, jj, kk) = true;
        % 
        %     hg2 = patch(isosurface(Xm_, Ym_, Zm_, G | single_voxel, 0.5));
        %     set(hg2,'EdgeColor','none','FaceColor',[0.1 0.6 1.0],'FaceAlpha',0.25);
        %     V0 = hg2.Vertices - 1; 
        %     Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
        %     hg2.Vertices = Vras(:,1:3);
        % end
        % % draw new tube (if available)
        % if isfinite(t0_new) && isfinite(t1_new) && isfinite(r_new)
        %     local_draw_tube(ctr_new, dstar_new, t0_new, t1_new, r_new);
        % end
        % scatter3(cand_ras(1),cand_ras(2),cand_ras(3),36,'r','filled');
        % axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
        % xlabel X; ylabel Y; zlabel Z;
        % title('Blob + candidate — refit tube');
        % % (4) Text panel: fill metrics
        
    
        % decide
        if (dFill < -fill_drop_tol) || (score < ScoreMin) || (fill_new < FillMin)
            decision = "reject";
            rejected(ii,jj,kk) = true;
            reason = " ";
            if (dFill < -fill_drop_tol)
                reason = [reason, 'too much fill reduction'];
            end
            if (score < ScoreMin)
                reason = [reason, 'low composite score'];
            end
            if (fill_new < FillMin)
                reason = [reason, 'too little fill ratio'];
            end
        else
            decision = "accept";
            reason = " ";
            G(ii,jj,kk) = true;   % accept
            grown = grown + 1;
            tube_change = true;
            nn_change = true;
            if grown >= max_voxels, stop_cause='cap_reached'; break; end
        end

        % continuous plot

        % nexttile; axis off;
        % txt = {
        %     sprintf('Decision: %s', decision)
        %     sprintf('reason: %s', reason)
        %     sprintf('cur fill = %.3f', cur_fill)
        %     sprintf('new fill = %.3f', fill_new)
        %     sprintf('Δfill = %.3f', dFill)
        %     sprintf('score = %.3f  (α=%.2f, β=%.2f)', score, W(1), W(2))
        %     sprintf('d_{min}=%.2f mm, d_{ave}=%.2f mm, pen=%.3f', dmin, dave, dist_pen)
        %     sprintf('tube now: r=%.2f mm, L=%.2f mm', r_mm, max(0,t1-t0))
        % };
        % text(0.01, 0.98, txt, 'VerticalAlignment','top','FontName','Monospaced');
        % 
        % sgtitle(sprintf('Grow step %d', grown));

        [logs, log_n] = append_log(logs, log_n, struct( ...
            'step',     grown, ...
            'ijk',      [ii jj kk], ...
            'ras',      cand_ras, ...
            'decision', decision, ...
            'reason', reason, ...
            'dFill',    dFill, ...
            'fill_new', fill_new, ...
            'cur_fill', cur_fill, ...
            'score',    score, ...
            'dist_pen', dist_pen, ...
            'dmin',     dmin, ...
            'dave',     dave, ...
            't0',       t0, ...
            't1',       t1, ...
            'r_mm',     r_mm, ...
            'ctr',      ctr, ...
            'dstar',    dstar, ...
            'note',     sprintf('KNearest=%d; AngleDeg=%.1f', KNearest, AngleDeg) ...
        ));

    
    end

    % ---------- pack info ----------
    info = struct( ...
        'grown_count', grown, ...
        'stop_cause',  stop_cause, ...
        'reason_counts', rc, ...
        'seed_ijk', seed_ijk, ...
        'AngleDeg', AngleDeg, ...
        'DirGrid', DirGrid, ...
        'Weights', W, ...
        'ScoreMin', ScoreMin, ...
        'KNearest', KNearest, ...
        'cur_fill', cur_fill);

    % figure('Color','w'); hold on;
    % % Blob surface (G) in RAS
    % if any(G(:))
    %     [m,n,p] = size(G);
    %     [Xm,Ym,Zm] = meshgrid(1:n,1:m,1:p);
    %     hg = patch(isosurface(Xm,Ym,Zm,G,0.5));
    %     set(hg,'EdgeColor','none','FaceColor',[1 0.7 0],'FaceAlpha',0.3);
    %     V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
    %     hg.Vertices = Vras(:,1:3);
    % end
    % 
    % % Planned line (span length based on voxel bbox extent)
    % U0 = U_in(:).'/max(norm(U_in),eps);
    % % quick length: use G extent if available; else a default span
    % if any(G(:))
    %     [ii,jj,kk] = ind2sub(size(G), find(G));
    %     PG = [ squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*1))), ...
    %            squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*2))), ...
    %            squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*3))) ];
    %     tG = (PG - start_traj) * U0(:);
    %     L0 = [min(tG) max(tG)] + [-TPadMM TPadMM];
    % else
    %     L0 = [-10 10]; % mm fallback
    % end
    % P0 = start_traj + L0(1)*U0; P1 = start_traj + L0(2)*U0;
    % plot3([P0(1) P1(1)],[P0(2) P1(2)],[P0(3) P1(3)],'-','LineWidth',2,'Color',[0 0.45 0.74]);
    % 
    % ctr = blob_centroid_ras(G, Xras); %
    % 
    % [dstar, cur_fill, t0, t1, r_mm] = fit_axis_by_fill( ...
    %     U, ctr, TPadMM, G, Xras, vx, ...
    %     'AngleDeg', AngleDeg, 'RingK', 12, 'IncludeCenter', true, ...
    %     'RGuardMM', 0.2, 'TGuardMM', 0.5, 'RobustPct', 98);
    % 
    % [t0, t1, r_mm, stats] = enclosing_tube_params(G, Xras, ctr, dstar);
    % 
    % scatter3(ctr(1),ctr(2),ctr(3),80,[0.85 0.33 0.10],'filled','MarkerEdgeColor','k');
    % 
    % % longitudinal sampling (denser if longer tube)
    % L  = max(t1 - t0, 1e-3);
    % nTau = max(24, ceil(L/3));     % ~3 mm spacing along axis
    % nTh  = 48;                      % circular resolution
    % tmp = [1;0;0];
    % if abs(dot(tmp,dstar)) > 0.99
    %     tmp = [0;1;0];
    % end
    % n1 = cross(dstar, tmp); n1 = n1 / norm(n1);
    % n2 = cross(dstar, n1);  n2 = n2 / norm(n2);
    % tau   = linspace(t0, t1, nTau);
    % theta = linspace(0, 2*pi, nTh);
    % [TH, TAU] = meshgrid(theta, tau);
    % 
    % % parametric tube surface in RAS
    % Xt = ctr(1) + TAU.*dstar(1) + r_mm*(cos(TH)*n1(1) + sin(TH)*n2(1));
    % Yt = ctr(2) + TAU.*dstar(2) + r_mm*(cos(TH)*n1(2) + sin(TH)*n2(2));
    % Zt = ctr(3) + TAU.*dstar(3) + r_mm*(cos(TH)*n1(3) + sin(TH)*n2(3));
    % 
    % hs_tube = surf(Xt, Yt, Zt, ...
    %     'FaceAlpha', 0.10, 'EdgeColor', 'none', 'FaceColor', [0 0.4 1]);
    % 
    % % (optional) end-caps to make the tube extent obvious
    % capTh = linspace(0, 2*pi, 64);
    % C0 = ctr + t0*dstar + r_mm*(cos(capTh)'*n1 + sin(capTh)'*n2);
    % C1 = ctr + t1*dstar + r_mm*(cos(capTh)'*n1 + sin(capTh)'*n2);
    % patch('XData',C0(:,1),'YData',C0(:,2),'ZData',C0(:,3), ...
    %       'FaceColor',[0 0.4 1],'FaceAlpha',0.10,'EdgeColor','none');
    % patch('XData',C1(:,1),'YData',C1(:,2),'ZData',C1(:,3), ...
    %       'FaceColor',[0 0.4 1],'FaceAlpha',0.10,'EdgeColor','none');
    % 
    % axis equal vis3d; grid on; box on; camlight headlight; lighting gouraud;
    % xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
    % title('Planned trajectory (blue), fitted tube (blue translucent), grown blob (orange)');
    % view(135,20);
    % 
    % figure('Color','w'); hold on;
    % 
    % hg = patch(isosurface(Xm,Ym,Zm,metal_roi_lo,0.5));
    % set(hg,'EdgeColor','none','FaceColor',[0 0 1],'FaceAlpha',0.3);
    % V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
    % hg.Vertices = Vras(:,1:3);
    % axis equal vis3d; grid on; box on; camlight headlight; lighting gouraud;
    % xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
    % title('Planned trajectory (blue), fitted tube (blue translucent), grown blob (orange)');
    % view(135,20);

end