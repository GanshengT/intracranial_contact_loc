function out = fit_shank_line_from_blob(metal_roi_unique, BW_skull, Avox2ras0, start_traj, end_traj, varargin)
    % Fit a shank line to a selected metal∩ROI blob (CT grid), clipped by brain mask.
    % Returns proximal_entry_ras (at brain boundary), distal_ras (from metal tail), dir_unit, etc.
    
    % ---- Options ----
    p = inputParser;
    p.addParameter('StepMM', 0.2);
    p.addParameter('TubeRadiusMM', 0.8);
    p.addParameter('Refine', false);              % default off; distal-from-cloud is robust
    p.addParameter('MaxDistalGrow', 5);
    p.addParameter('BrainMaskCT', []);            % logical, CT grid
    p.addParameter('BrainMaskMGZ', []);           % struct from MRIread or numeric vol
    p.addParameter('BrainMaskVox2RAS', []);       % MRI vox->RAS
    p.addParameter('VisualizeBands', false);
    p.addParameter('VisualizeFinal', false);
    p.addParameter('VisualizeTube_rmm', false);
    p.addParameter('VisualizeTube_rmm_975', false);
    p.addParameter('subj_id', []);
    p.addParameter('report_dir', []);
    p.addParameter('traj_id', []);
    p.addParameter('curved', false);
    p.parse(varargin{:});
    opt = p.Results;
    
    % ---- 0) Brain mask in CT space ----
    BW_brain_CT = [];
    if ~isempty(opt.BrainMaskCT)
        assert(isequal(size(opt.BrainMaskCT), size(metal_roi_unique)), 'BrainMaskCT size mismatch.');
        BW_brain_CT = logical(opt.BrainMaskCT);
    elseif ~isempty(opt.BrainMaskMGZ) && ~isempty(opt.BrainMaskVox2RAS)
        if isstruct(opt.BrainMaskMGZ) && isfield(opt.BrainMaskMGZ,'vol')
            BW_brain_MRI = opt.BrainMaskMGZ.vol > 0;
        else
            BW_brain_MRI = logical(opt.BrainMaskMGZ);
        end
        A_mri = opt.BrainMaskVox2RAS;    % MRI vox(0-based)->RAS
        BW_brain_CT = resample_mask_nn(BW_brain_MRI, A_mri, size(metal_roi_unique), Avox2ras0);
    end
    if isempty(BW_brain_CT)
        % fall back to skull-interior (less robust than true brain mask)
        interior = imfill(BW_skull,'holes') & ~BW_skull;
        BW_brain_CT = interior;
    end
    % diagnostic: check brain mask from frresurfer
    % results: passed, isosurface as intended
    % [Xm, Ym, Zm] = meshgrid(1:size(BW_brain_CT,1), ...
    %                   1:size(BW_brain_CT,2), ...
    %                   1:size(BW_brain_CT,3));
    % figure('Color','w'); hold on;
    % 
    % % --- Plot skull (outer) surface ---
    % if any(BW_skull(:))
    %     skull_iso = isosurface(Xm, Ym, Zm, BW_skull, 0.5);
    %     V0 = skull_iso.vertices - 1;  % convert to 0-based before RAS transform
    %     Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.';
    %     skull_iso.vertices = Vras(:,1:3);
    %     h_skull = patch(skull_iso);
    %     set(h_skull, 'EdgeColor','none', 'FaceColor', [0.9 0.7 0.3], 'FaceAlpha', 0.25);
    % end
    % 
    % % --- Plot brain mask surface ---
    % if any(BW_brain_CT(:))
    %     brain_iso = isosurface(Xm, Ym, Zm, BW_brain_CT, 0.5);
    %     V0 = brain_iso.vertices - 1;  % convert to 0-based
    %     Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.';
    %     brain_iso.vertices = Vras(:,1:3);
    %     h_brain = patch(brain_iso);
    %     set(h_brain, 'EdgeColor','none', 'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.4);
    % end
    % axis equal vis3d
    % camlight headlight; lighting gouraud
    % xlabel('X_RAS (mm)'); ylabel('Y_RAS (mm)'); zlabel('Z_RAS (mm)');
    % legend([h_skull h_brain], {'Skull', 'Brain Mask'}, 'Location','best');
    % title('Skull and Brain Mask in RAS space');
    
    % ---- 1) Split metal by brain mask & pick inside for fitting ----
    % note that BW_brain_CT is unlikely have holes (because we filled it)
    metal_in  = metal_roi_unique & BW_brain_CT;
    metal_out = metal_roi_unique & ~BW_brain_CT;
    inside_ct  = nnz(metal_in);
    outside_ct = nnz(metal_out);
    assert(any(metal_in(:)), 'No metal voxels remain inside brain mask/skull interior.');
    
    % ---- 2) Collect inside-metal points in RAS ----
    [idx_i, idx_j, idx_k] = ind2sub(size(metal_in), find(metal_in));
    % IMPORTANT: map [i j k] -> [x y z] = [j i k] before Avox2ras
    % THIS function ask for Avox2ras0
    P = ijk1_to_ras(Avox2ras0, [idx_j idx_i idx_k]);   % Nx3 RAS
    
    % ---- 3) Robust axis fit (no planned line) ----
    if ~opt.curved
        [rfit, axis_stats] = fit_axis_cylinder_robust(P, ...
            'MaxIter', 20, 'Tol', 1e-5, 'RadiusGrid', linspace(0.3, 2.0, 9), 'subj_id', opt.subj_id, ...
            'report_dir', opt.report_dir, 'traj_id', opt.traj_id);
            p0  = rfit.p0;          % RAS
            e1  = rfit.u;  e1 = e1 ./ max(norm(e1), eps);
            rmm = rfit.r; % best‑fit cylinder radius in mm 
    else
        res = fit_axis_piecewise_knee(P, ...
            'DetectKnee', true, 'Refine', true, ...
            'Report', struct('subj_id', opt.subj_id, 'report_dir', opt.report_dir, 'traj_id', opt.traj_id), ...
            'Visualize', true);
    end


    
    % ---- 4) Proximal from brain boundary; orient e1 inward ----
    if opt.curved
        % --- Curved mode ---
        % Determine which parts of the fitted spline are inside the brain mask
        C_curve = res.centerline;     % M×3 centerline points
        P_all   = P;                  % voxel cloud
        M       = size(C_curve,1);
    
        % --- voxel-to-centerline distances ---
        % For each voxel, find the nearest point on the curve
        % (use vectorized nearest-neighbor approach)
        D = pdist2(P_all, C_curve);      % N×M distances
        [dmin, idx_min] = min(D, [], 2); % nearest centerline index per voxel
    
        % --- tube radius ---
        % define inliers as voxels within some quantile distance
        r_best = prctile(dmin, 90);              % 90th percentile (analogous to Tukey in straight fit)
        in_best = dmin <= r_best;                % voxel mask inside curved tube
        r_median = median(dmin(in_best));
    
        % (optional) radius profile along the curve (for visualization)
        r_of_t = zeros(M,1);
        for j = 1:M
            r_of_t(j) = median(dmin(idx_min==j));  % median radius at each centerline node
        end
    
        % --- (3) Begin / end and span along curve ---
        % Compute arc-length cumulative distance
        seglen = sqrt(sum(diff(C_curve,1,1).^2, 2));
        s_curve = [0; cumsum(seglen)];
        total_span = s_curve(end);
    
        % Assign start and end points
        proximal_entry_ras = C_curve(1,:);
        distal_ras         = C_curve(end,:);
        t_span_mm          = total_span;
    
        % --- (4) Optional: orientation correction based on brain center ---
        if exist('BW_brain_CT','var') && ~isempty(BW_brain_CT)
            C_brain = brain_centroid_ras(BW_brain_CT, Avox2ras0);
            M_mid = 0.5*(proximal_entry_ras + distal_ras);
            dir_vec = distal_ras - proximal_entry_ras;
            if dot(dir_vec, C_brain - M_mid) < 0
                % flip the direction if it points outward
                C_curve = flipud(C_curve);
                proximal_entry_ras = C_curve(1,:);
                distal_ras         = C_curve(end,:);
                s_curve = s_curve(end) - s_curve;  % reverse param
            end
        end
    
        % --- (5) Pack results into rfit-like struct ---
        rfit = struct( ...
            'p0', C_curve(1,:), ...
            'u', (C_curve(end,:) - C_curve(1,:)) ./ max(norm(C_curve(end,:) - C_curve(1,:)), eps), ...
            'r', r_best, ...
            'r_median', r_median, ...
            'r_of_t', r_of_t, ...
            'centerline', C_curve, ...
            'inliers', in_best, ...
            'entry', proximal_entry_ras, ...
            'exit', distal_ras, ...
            'span', t_span_mm ...
        );
        % figure('Color','w','Position',[100 100 950 700]); hold on; axis equal vis3d;
        % scatter3(P_all(:,1), P_all(:,2), P_all(:,3), 8, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.25);
        % plot3(C_curve(:,1), C_curve(:,2), C_curve(:,3), 'k-', 'LineWidth', 2.0);
        % scatter3(proximal_entry_ras(1), proximal_entry_ras(2), proximal_entry_ras(3), 40, 'g', 'filled');
        % scatter3(distal_ras(1), distal_ras(2), distal_ras(3), 40, 'r', 'filled');
        % title(sprintf('Curved axis fit: radius=%.2f mm, span=%.1f mm', r_best, t_span_mm));
        % legend({'Voxels','Centerline','Entry','Exit'}, 'Location','best');
        % xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)'); grid on; view(135,25);
        % filename = [opt.subj_id, '_curved_estimation_entry_and_end_', opt.traj_id];
        % full_path = fullfile(opt.report_dir, [filename, '.pdf']);
        % saveas(gcf, full_path, 'pdf');
        rmm = r_best;
        model_chosen = 'polyline';
        centerline_pts = C_curve;
    
    else
        % --- Straight-line mode (original code) ---
        step_mm = max(0.2, opt.StepMM);
        pad_mm  = 5;
        t_all = (P - p0) * e1.';
        t_span  = [min(t_all)-pad_mm, max(t_all)+pad_mm];
    
        [B1, B2, t1, t2] = brain_boundary_hits(BW_brain_CT, Avox2ras0, p0, e1, t_span, step_mm);
        assert(all(isfinite([t1 t2])), 'Boundary hits not found; check mask/affine/axis.');
    
        C = brain_centroid_ras(BW_brain_CT, Avox2ras0);
        M = 0.5*(B1 + B2);
        if dot(e1, C - M) < 0
            e1 = -e1; [B1,B2] = deal(B2,B1); [t1,t2] = deal(t2,t1);
        end
        proximal_entry_ras = B1;
        t2_rel = (B2 - proximal_entry_ras) * e1.';  % brain limit along +e1
    
        % Distal point along the fitted axis within the brain
        distal_ras = distal_from_pointcloud_robust( ...
            P, proximal_entry_ras, e1, ...
            'BinMM', 0.25, 'SmoothWin', 9, 'TailPct', 98.5, ...
            'ClipTMax', t2_rel, 'RadQuantile', 0.90);
    
        % Optional refinement (rarely needed)
        if opt.Refine
            [distal_ras_tmp, score_tmp] = refine_distal_tube(P, proximal_entry_ras, e1, distal_ras, opt.TubeRadiusMM, opt.MaxDistalGrow);
            if ((distal_ras_tmp - proximal_entry_ras) * e1.') <= (t2_rel + 1e-6)
                distal_ras = distal_ras_tmp;
            end
            score = score_tmp;
        else
            t_from_prox = (P - proximal_entry_ras) * e1.';
            score = mean(t_from_prox >= 0 & t_from_prox <= ((distal_ras - proximal_entry_ras) * e1.'));
        end
    
        dir_unit = distal_ras - proximal_entry_ras;
        dir_unit = dir_unit ./ max(norm(dir_unit), eps);
        t_span_mm = norm(distal_ras - proximal_entry_ras);
    end


    
    % possibility for adding a third pt and use spline intepolation,
    % see if spline interpolation is better than linear fitting (could
    % determine by baysien, like information criterion, is it worthy to add one parameter)

    E = proximal_entry_ras; 
    D = distal_ras;
    N = size(P,1);
    % Choose Rmax sensibly; 1.5–2.0× your tube radius works well
    Rmax_mm = min(2 * 0.8, rmm);
    if opt.curved
        % --- Curved: compare Bézier vs. actual polyline ---
        % Polyline distance & AUC
        d_poly  = sqrt(dist2_point_polyline(P, centerline_pts));
        RSS_poly = robust_rss(d_poly,5);
        AUC_poly = tube_auc_score_general(P, 'polyline', struct('C', centerline_pts), Rmax_mm, 5);
        BIC_poly = N * log(max(RSS_poly/N,eps)) + 3*log(N);  % small penalty
        DeltaBIC = 0; BIC_line = NaN; BIC_bez = NaN; AUC_line_curve = NaN; AUC_bez_curve = NaN;
        winner = "Polyline (curved fit)";
        use_bezier = false;
    else
        % --- Straight: compare line vs Bézier (your existing code block) ---
        [Ctrl_bez, RSS_bez] = fit_bezier_control_point(P, E, D, ...
            'SamplesM', 201, 'TrimPct', 5, 'Lambda', 0.0);
        d_line = sqrt(dist2_point_line_segment(P, E, D));

        RSS_line = robust_rss(d_line, 5);
        d_bez = sqrt(dist2_point_bezier(P, E, Ctrl_bez, D, 201));
        RSS_bez = robust_rss(d_bez, 5);
        BIC_line = N*log(max(RSS_line/N,eps));
        BIC_bez  = N*log(max(RSS_bez/N,eps))+3*log(N);
        DeltaBIC = BIC_line - BIC_bez;
        AUC_line_curve = tube_auc_score_general(P,'line',struct('E',E,'D',D),Rmax_mm,5);
        AUC_bez_curve  = tube_auc_score_general(P,'bezier',struct('E',E,'C',Ctrl_bez,'D',D,'SamplesM',201),Rmax_mm,5);
        winner = ternary(DeltaBIC>0,'Bézier','Line');
        use_bezier = (DeltaBIC > 10);
    end

  
    if opt.VisualizeTube_rmm
        nTheta = 40; nSeg = 40; % mesh resolution
        figure('Color','w'); hold on;
        % [Xm, Ym, Zm] = meshgrid(1:size(BW_brain_CT,1), ...
    %                   1:size(BW_brain_CT,2), ...
    %                   1:size(BW_brain_CT,3));
    
        % (1) Visualize base metal structures
        [Xm, Ym, Zm] = meshgrid(1:size(metal_in,2),1:size(metal_in,1),1:size(metal_in,3));
        toRAS = @(V1) ([V1-1, ones(size(V1,1),1)] * Avox2ras0.');
    
        if any(metal_in(:))
            hi = patch(isosurface(Xm,Ym,Zm,metal_in,0.5));
            set(hi,'EdgeColor','none','FaceAlpha',0.65,'FaceColor',[0.90 0.20 0.20]);
            V0 = hi.Vertices; Vras = toRAS(V0); hi.Vertices = Vras(:,1:3);
        end
        if any(metal_out(:))
            ho = patch(isosurface(Xm,Ym,Zm,metal_out,0.5));
            set(ho,'EdgeColor','none','FaceAlpha',0.60,'FaceColor',[0.10 0.45 0.85]);
            V0 = ho.Vertices; Vras = toRAS(V0); ho.Vertices = Vras(:,1:3);
        end
    
        % (2) Tube visualization depending on model
        if opt.curved
            % ---- Curved centerline (polyline) ----
            Cseg = centerline_pts;
            % Draw continuous small tube segments along the curve
            for k = 1:size(Cseg,1)-1
                [Vb, Fb] = tube_mesh_segment(Cseg(k,:), Cseg(k+1,:), rmm, nTheta, 1);
                if ~isempty(Vb)
                    patch('Vertices', Vb, 'Faces', Fb, ...
                          'FaceColor', [0.10 0.85 0.45], ...
                          'EdgeColor', 'none', 'FaceAlpha', 0.15);
                end
            end
            % Plot the centerline itself
            plot3(Cseg(:,1), Cseg(:,2), Cseg(:,3), 'k-', 'LineWidth', 2);
            scatter3(Cseg(1,1),Cseg(1,2),Cseg(1,3),40,'g','filled');
            scatter3(Cseg(end,1),Cseg(end,2),Cseg(end,3),40,'r','filled');
            title(sprintf('Curved lead fit: tube r=%.2f mm', rmm));
    
        elseif exist('Ctrl_bez','var') && ~isempty(Ctrl_bez)
            % ---- Bézier centerline ----
            [~, Bcen] = bezier_samples(E, Ctrl_bez, D, 60);
            for k = 1:size(Bcen,1)-1
                [Vb, Fb] = tube_mesh_segment(Bcen(k,:), Bcen(k+1,:), rmm, nTheta, 1);
                if ~isempty(Vb)
                    patch('Vertices', Vb, 'Faces', Fb, ...
                          'FaceColor', [0.10 0.85 0.45], ...
                          'EdgeColor', 'none', 'FaceAlpha', 0.15);
                end
            end
            plot3(Bcen(:,1),Bcen(:,2),Bcen(:,3),'g-','LineWidth',2);
            scatter3(E(1),E(2),E(3),40,'g','filled');
            scatter3(D(1),D(2),D(3),40,'r','filled');
            title(sprintf('Bézier lead fit: tube r=%.2f mm', rmm));
    
        else
            % ---- Straight line segment ----
            [Vt, Ft] = tube_mesh_segment(E, D, rmm, nTheta, nSeg);
            if ~isempty(Vt)
                patch('Vertices', Vt, 'Faces', Ft, ...
                      'FaceColor', [0.10 0.45 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.18);
            end
            plot3([E(1) D(1)], [E(2) D(2)], [E(3) D(3)], 'k-', 'LineWidth', 2);
            scatter3(E(1),E(2),E(3),40,'g','filled');
            scatter3(D(1),D(2),D(3),40,'r','filled');
            title(sprintf('Straight line fit: tube r=%.2f mm', rmm));
        end
    
        % (3) In-tube point visualization (works for all models)
        if opt.curved
            d2_vis = dist2_point_polyline(P, centerline_pts);
        elseif exist('Ctrl_bez','var') && ~isempty(Ctrl_bez)
            d2_vis = dist2_point_bezier(P, E, Ctrl_bez, D, 201);
        else
            d2_vis = dist2_point_line_segment(P, E, D);
        end
        inTube = sqrt(d2_vis) <= rmm;
    
        % Subsample for clarity
        nPlot = min(5000, size(P,1));
        idx = randperm(size(P,1), nPlot);
        idx_in  = idx(inTube(idx));
        idx_out = idx(~inTube(idx));
        scatter3(P(idx_out,1), P(idx_out,2), P(idx_out,3), ...
            6, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.25);
        scatter3(P(idx_in,1), P(idx_in,2), P(idx_in,3), ...
            6, [0.1 0.7 0.1], 'filled', 'MarkerFaceAlpha', 0.45);
    
        % (4) Global aesthetics
        axis equal vis3d; camlight headlight; lighting gouraud;
        xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
        view(135,20);
        filename = [opt.subj_id, '_modeling_lead_' , opt.traj_id, '_tube_visualization'];
        full_path = fullfile(opt.report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');
    end

    % diagnostic for AUC calculation
    trimPct = 5;
    R = linspace(0, Rmax_mm, 256);
    
    % Prepare containers
    AUC_line_curve = NaN; AUC_bez_curve = NaN; AUC_poly_curve = NaN;
    F_line = []; F_bez = []; F_poly = [];
    
    if opt.curved
        % --- Polyline / curved centerline ---
        d_poly = sqrt(dist2_point_polyline(P, centerline_pts));
        d_poly = sort(d_poly);
        n = numel(d_poly);
        k0 = floor(n*trimPct/100)+1; k1 = n - floor(n*trimPct/100);
        d_poly = d_poly(max(1,k0):max(k0,k1));
        F_poly = arrayfun(@(r) mean(d_poly <= r), R);
        AUC_poly_curve = trapz(R, F_poly) / Rmax_mm;
    
        % --- Plot ---
        % fig_auc = figure('Color','w'); hold on;
        % plot(R, F_poly, 'b-', 'LineWidth', 2);
        % grid on; box on;
        % xlabel('Radius r (mm)'); ylabel('Fraction within r (CDF)');
        % title('Tube AUC curve (curved polyline fit)');
        % legend({'Polyline'}, 'Location','southeast');
        % text(0.6*Rmax_mm, 0.15, sprintf('AUC_{poly}=%.3f  BIC=%.1f', ...
        %     AUC_poly_curve, BIC_poly), 'Color','b');
    
        % Model tag
        model_chosen = 'polyline';
        use_bezier = false;
        t_samp = linspace(0,1,200).';
        % Normalize centerline sampling along arc length if needed
        centerline_pts = resample_centerline_by_arclength(centerline_pts, numel(t_samp));
    
    else
        % --- Line ---
        d_line = sqrt(dist2_point_line_segment(P, E, D));
        d_line = sort(d_line);
        n = numel(d_line);
        k0 = floor(n*trimPct/100)+1; k1 = n - floor(n*trimPct/100);
        d_line = d_line(max(1,k0):max(k0,k1));
        F_line = arrayfun(@(r) mean(d_line <= r), R);
        AUC_line_curve = trapz(R, F_line) / Rmax_mm;
    
        % --- Bézier ---
        d_bez = sqrt(dist2_point_bezier(P, E, Ctrl_bez, D, 201));
        d_bez = sort(d_bez);
        n2 = numel(d_bez);
        k0b = floor(n2*trimPct/100)+1; k1b = n2 - floor(n2*trimPct/100);
        d_bez = d_bez(max(1,k0b):max(k0b,k1b));
        F_bez = arrayfun(@(r) mean(d_bez <= r), R);
        AUC_bez_curve = trapz(R, F_bez) / Rmax_mm;
    
        % --- Plot comparison ---
        % fig_auc = figure('Color','w'); hold on;
        % plot(R, F_line, 'k-', 'LineWidth', 2);
        % plot(R, F_bez,  'g-', 'LineWidth', 2);
        % grid on; box on;
        % xlabel('Radius r (mm)'); ylabel('Fraction within r (CDF)');
        % title('Tube AUC curves (trimmed) — Line vs Bézier');
        % 
        [bestBIC, bestIdx] = min([BIC_line, BIC_bez]);
        dBIC_line = BIC_line - bestBIC;
        dBIC_bez  = BIC_bez  - bestBIC;
        winner = ["Line","Bézier"];
        winner = winner(bestIdx);
    
        % comment for masking the plots
        % leg = legend({'Line','Bézier'}, 'Location','southeast');
        % set(leg, 'AutoUpdate','off');
        % txty = 0.15;
        % text(0.60*Rmax_mm, txty+0.16, sprintf('Winner: %s (lower BIC)', winner), ...
        %     'FontWeight','bold', 'Color','k');
        % text(0.60*Rmax_mm, txty+0.10, sprintf('Line   AUC=%.3f  BIC=%.1f  ΔBIC=%.1f ', ...
        %     AUC_line_curve, BIC_line, dBIC_line), 'Color','k');
        % text(0.60*Rmax_mm, txty+0.04, sprintf('Bézier AUC=%.3f  BIC=%.1f  ΔBIC=%.1f ', ...
        %     AUC_bez_curve,  BIC_bez,  dBIC_bez),  'Color','g');
    
        % Optional: star the winner curve
        % if bestIdx == 1
        %     plot(NaN,NaN,'k*','MarkerSize',10,'DisplayName','Winner');
        % else
        %     plot(NaN,NaN,'g*','MarkerSize',10,'DisplayName','Winner');
        % end
        % legend('show');
    
        use_bezier = (DeltaBIC > 10);
        t_samp = linspace(0,1,200).';
        if use_bezier
            [~, centerline_pts] = bezier_samples(E, Ctrl_bez, D, numel(t_samp));
            model_chosen = 'bezier';
        else
            centerline_pts = E + t_samp.*(D - E);
            model_chosen = 'line';
        end
    end
    
    % --- Save diagnostic figure --- optional, commented
    % filename = [opt.subj_id, '_tube_auc_diagnostic_', opt.traj_id];
    % full_path = fullfile(opt.report_dir, [filename, '.pdf']);
    % saveas(gcf, full_path, 'pdf');
    
   
    % Use inliers if available; else all P
    if isfield(rfit,'inliers') && ~isempty(rfit.inliers)
        Pin = P(rfit.inliers, :);
    else
        Pin = P;
    end
    
    % Radial distances to the chosen centerline
    % Axial trimming to avoid fringe artifacts (same trimPct you use elsewhere)
    trimPct = 5;

    % --- Choose proper distance & axis projection functions ---
    if opt.curved
        % Distances from polyline centerline
        d_rad = sqrt(dist2_point_polyline(Pin, centerline_pts));   % mm
        
        % Define pseudo-axial coordinate (arc length projection)
        seglen = sqrt(sum(diff(centerline_pts,1,1).^2,2));
        s_curve = [0; cumsum(seglen)];
        % Project each point to nearest curve location
        [~, idx_near] = min(pdist2(Pin, centerline_pts), [], 2);
        t_ax = s_curve(idx_near);
    else
        % Straight/Bézier
        d_rad = radial_dist_to_centerline(Pin, E, D, model_chosen, Ctrl_bez);  % mm
        t_ax  = project_to_axis(Pin, E, D);
    end

    lo = prctile(t_ax, trimPct);
    hi = prctile(t_ax, 100 - trimPct);
    keep = (t_ax >= lo) & (t_ax <= hi);
    d_use = d_rad(keep);
    t_use = t_ax(keep);
    
    % Global quantiles
    r_q = prctile(d_use, [50 80 90 95 97.5]);   % mm
    
    % Optional: per-axial-bin conservative 95th profile (catches local widenings)
    axBinMM = 2.0;
    edges = min(t_use):axBinMM:max(t_use);
    r95_bins = nan(max(numel(edges)-1,0),1);
    for b = 1:numel(edges)-1
        idxb = (t_use >= edges(b)) & (t_use < edges(b+1));
        if nnz(idxb) > 10
            r95_bins(b) = prctile(d_use(idxb), 95);
        end
    end
    r95_profile_max = max(r95_bins, [], 'omitnan');   % conservative
    if isfinite(r95_profile_max)
        r_q(4) = r95_profile_max;  % replace global 95th with per-bin max 95th
    end
    
    % Name them for clarity
    rmm_50  = r_q(1);
    rmm_80  = r_q(2);
    rmm_90  = r_q(3);
    rmm_95  = r_q(4);
    rmm_975 = r_q(5);

    if opt.VisualizeTube_rmm_975
        nTheta = 40; nSeg = 40; % mesh resolution
        figure('Color','w'); hold on;
    
        if opt.curved
            % ===== Polyline curved case =====
            % Draw tube segments along the fitted curved centerline
            Cseg = centerline_pts;
            for k = 1:size(Cseg,1)-1
                [Vb, Fb] = tube_mesh_segment(Cseg(k,:), Cseg(k+1,:), rmm_975, nTheta, 1);
                if ~isempty(Vb)
                    patch('Vertices', Vb, 'Faces', Fb, ...
                          'FaceColor',[0.10 0.85 0.45], ...
                          'EdgeColor','none','FaceAlpha',0.18);
                end
            end
            % Centerline and endpoints
            plot3(Cseg(:,1),Cseg(:,2),Cseg(:,3),'k-','LineWidth',2);
            scatter3(Cseg(1,1),Cseg(1,2),Cseg(1,3),40,'g','filled');
            scatter3(Cseg(end,1),Cseg(end,2),Cseg(end,3),40,'r','filled');
            title(sprintf('Curved polyline fit: r=%.2f mm', rmm_975));
    
            % (Optional) visualize in-tube points
            d2_vis = dist2_point_polyline(P, Cseg);
            inTube = sqrt(d2_vis) <= rmm_975;
            nPlot = min(5000, size(P,1));
            idx = randperm(size(P,1), nPlot);
            idx_in  = idx(inTube(idx));
            idx_out = idx(~inTube(idx));
            scatter3(P(idx_out,1), P(idx_out,2), P(idx_out,3), ...
                6, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.25);
            scatter3(P(idx_in,1), P(idx_in,2), P(idx_in,3), ...
                6, [0.1 0.7 0.1], 'filled', 'MarkerFaceAlpha', 0.45);
    
        else
            % ===== Original (line + Bézier) code block unchanged =====
            [Eseg, Dseg] = deal(proximal_entry_ras, distal_ras);
            [Vt, Ft] = tube_mesh_segment(Eseg, Dseg, rmm_975, nTheta, nSeg);
            if ~isempty(Vt)
                htube = patch('Vertices', Vt, 'Faces', Ft, ...
                    'FaceColor', [0.10 0.45 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.18);
            end
            
            % (Optional) visualize which metal-in points fall inside the tube
            d2_line_vis = dist2_point_line_segment(P, Eseg, Dseg);
            inTube = sqrt(d2_line_vis) <= rmm_975;
            % Subsample for speed/clarity
            nPlot = min(5000, size(P,1));
            idx = randperm(size(P,1), nPlot);
            idx_in  = idx(inTube(idx));
            idx_out = idx(~inTube(idx));
            
            scatter3(P(idx_out,1), P(idx_out,2), P(idx_out,3), ...
                6, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.25);
            
            scatter3(P(idx_in,1), P(idx_in,2), P(idx_in,3), ...
                6, [0.1 0.7 0.1], 'filled', 'MarkerFaceAlpha', 0.45);
            
            if any(metal_in(:))
                hi = patch(isosurface(Xm,Ym,Zm,metal_in,0.5));
                set(hi,'EdgeColor','none','FaceAlpha',0.65,'FaceColor',[0.90 0.20 0.20]);
                V0 = hi.Vertices; Vras = toRAS(V0); hi.Vertices = Vras(:,1:3);
            end
            if any(metal_out(:))
                ho = patch(isosurface(Xm,Ym,Zm,metal_out,0.5));
                set(ho,'EdgeColor','none','FaceAlpha',0.60,'FaceColor',[0.10 0.45 0.85]);
                V0 = ho.Vertices; Vras = toRAS(V0); ho.Vertices = Vras(:,1:3);
            end
    
            if exist('Ctrl_bez','var') && ~isempty(Ctrl_bez)
                [tB, Bcen] = bezier_samples(Eseg, Ctrl_bez, Dseg, 60);
                for k = 1:size(Bcen,1)-1
                    [Vb, Fb] = tube_mesh_segment(Bcen(k,:), Bcen(k+1,:), rmm_975, 24, 1);
                    if ~isempty(Vb)
                        patch('Vertices', Vb, 'Faces', Fb, ...
                              'FaceColor', [0.10 0.85 0.45], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
                    end
                end
                plot3(Bcen(:,1),Bcen(:,2),Bcen(:,3),'g-','LineWidth',2);
            end
        end
    
        % ===== Shared visual setup =====
        axis equal vis3d; camlight headlight; lighting gouraud;
        xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
        title(sprintf('Metal-in (red), Tube r=%.2f mm (blue), In-tube points (green)', rmm_975));
        view(135,20);
        filename = [opt.subj_id, '_modeling_lead_' , opt.traj_id, '_by_tube_fiting_rmm975'];
        full_path = fullfile(opt.report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');
    end

  
    if exist('B1','var') && exist('B2','var') && exist('t1','var') && exist('t2','var')
        boundary_hits = struct('B1',B1,'B2',B2,'t1',t1,'t2',t2);
    else
        boundary_hits = [];  % curved path or if boundary hits weren't computed
    end
    if ~exist('BIC_poly','var')
        BIC_poly = nan;
    end
    if ~exist('dir_unit', 'var')
        dir_unit = nan;
    end
    if ~exist('e1', 'var')
        e1 = nan;
    end
    if ~exist('score', 'var')
        score = nan;
    end
    if ~exist('axis_stats', 'var')
        axis_stats = nan;
    end
    if ~exist('Ctrl_bez', 'var')
        Ctrl_bez = [];
    end
    out = struct( ...
        'proximal_entry_ras', proximal_entry_ras, ...
        'distal_ras',         distal_ras, ...
        'dir_unit',           dir_unit, ...
        'pca_axis',           e1, ...
        't_span_mm',          t_span_mm, ...
        'score',              score, ...
        'clip_counts',        [inside_ct, outside_ct], ...
        'axis_fit',           rfit, ...
        'axis_stats',         axis_stats, ...
        'tube_radius_mm',     rmm, ...
        'boundary_hits',      boundary_hits, ...
        'tube_auc_line',      AUC_line_curve, ...
        'tube_auc_bezier',    AUC_bez_curve, ...
        'tube_auc_polyline',  AUC_poly_curve, ...
        'bic_line',           BIC_line, ...
        'bic_bezier',         BIC_bez, ...
        'bic_poly',           BIC_poly, ...
        'delta_bic',          DeltaBIC, ...
        'model',              model_chosen, ...
        'bezier_control',     Ctrl_bez, ...
        'centerline_t',       t_samp, ...
        'centerline_pts',     centerline_pts, ...
        'rmm_50',             rmm_50, ...
        'rmm_80',             rmm_80, ...
        'rmm_90',             rmm_90, ...
        'rmm_95',             rmm_95, ...
        'rmm_975',            rmm_975 );
end

function t = project_to_axis(P, E, D)
    u = (D - E); L = max(norm(u), eps); u = u / L;
    t = (P - E) * u.';  % mm along the axis from E
end

function d = radial_dist_to_centerline(P, E, D, model, Ctrl_bez)
    switch lower(model)
        case 'bezier'
            d = sqrt(dist2_point_bezier(P, E, Ctrl_bez, D, 201));  % mm
        otherwise
            d = sqrt(dist2_point_line_segment(P, E, D));           % mm
    end
end