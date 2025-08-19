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
    p.addParameter('VisualizeFinal', true);
    p.addParameter('VisualizeTube_rmm', true);
    p.addParameter('subj_id', []);
    p.addParameter('report_dir', []);
    p.addParameter('traj_id', []);
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
    
    % ---- 1) Split metal by brain mask & pick inside for fitting ----
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
    [rfit, axis_stats] = fit_axis_cylinder_robust(P, ...
        'MaxIter', 20, 'Tol', 1e-5, 'RadiusGrid', linspace(0.3, 2.0, 9));
    p0  = rfit.p0;          % RAS
    e1  = rfit.u;  e1 = e1 ./ max(norm(e1), eps);
    rmm = rfit.r;
    
    % Inlier-based axial bounds (for viz / span only)
    t_all = (P - p0) * e1.';
    if any(rfit.inliers)
        t_lo = prctile(t_all(rfit.inliers),  2.5);
        t_hi = prctile(t_all(rfit.inliers), 97.5);
    else
        t_lo = min(t_all); t_hi = max(t_all);
    end
    Plo = p0 + t_lo*e1;  Phi = p0 + t_hi*e1; 
    
    % ---- 4) Proximal from brain boundary; orient e1 inward ----
    step_mm = max(0.2, opt.StepMM);
    pad_mm  = 5;
    t_span  = [min(t_all)-pad_mm, max(t_all)+pad_mm];
    
    [B1, B2, t1, t2] = brain_boundary_hits(BW_brain_CT, Avox2ras0, p0, e1, t_span, step_mm);
    assert(all(isfinite([t1 t2])), 'Boundary hits not found; check mask/affine/axis.');
    
    C = brain_centroid_ras(BW_brain_CT, Avox2ras0);
    M = 0.5*(B1 + B2);
    if dot(e1, C - M) < 0
        e1 = -e1; [B1,B2] = deal(B2,B1); [t1,t2] = deal(t2,t1);
    end
    proximal_entry_ras = B1;
    t2_rel = (B2 - proximal_entry_ras) * e1.';    % brain limit along +e1
    
    % Distal from metal point-cloud tail (clipped to brain boundary) ----
    distal_ras = distal_from_pointcloud_robust( ...
        P, proximal_entry_ras, e1, ...
        'BinMM', 0.25, 'SmoothWin', 9, 'TailPct', 98.5, ...
        'ClipTMax', t2_rel, 'RadQuantile', 0.90);
    
    % Optional refine (kept off by default; ensure it stays inside brain if enabled)
    if opt.Refine
        [distal_ras_tmp, score_tmp] = refine_distal_tube(P, proximal_entry_ras, e1, distal_ras, opt.TubeRadiusMM, opt.MaxDistalGrow);
        % Clip back to boundary if refinement overshot
        if ((distal_ras_tmp - proximal_entry_ras) * e1.') <= (t2_rel + 1e-6)
            distal_ras = distal_ras_tmp;
        end
        score = score_tmp;
    else
        % Define a simple score without tube radius (fraction of points between ends)
        t_from_prox = (P - proximal_entry_ras) * e1.';
        score = mean(t_from_prox >= 0 & t_from_prox <= ( (distal_ras - proximal_entry_ras) * e1.' ));
    end
    
    % ---- 6) Direction & span (prox->distal) ----
    dir_unit = distal_ras - proximal_entry_ras;
    dir_unit = dir_unit ./ max(norm(dir_unit), eps);
    t_from_prox = (P - proximal_entry_ras) * e1.';
    t_span_mm   = max(t_from_prox) - max(0, min(t_from_prox));
    
    % --- viz (optional) ----
    if opt.VisualizeBands || opt.VisualizeFinal
        [m,n,p] = size(metal_roi_unique);
        [Xm, Ym, Zm] = meshgrid(1:n,1:m,1:p);
        toRAS = @(V1) ([V1-1, ones(size(V1,1),1)] * Avox2ras0.');
    
        if opt.VisualizeBands
            legend_handles = [];
            legend_labels = {};
            figure('Color','w'); hold on;
            if any(BW_skull(:))
                hs = patch(isosurface(Xm,Ym,Zm,BW_skull,0.5));
                set(hs,'EdgeColor','none','FaceColor',[0.85 0.80 0.70],'FaceAlpha',0.30);
                V0 = hs.Vertices; Vras = toRAS(V0); hs.Vertices = Vras(:,1:3);
                legend_handles(end+1) = hs;
                legend_labels{end+1} = 'Skull';
            end
            if any(metal_in(:))
                hi = patch(isosurface(Xm,Ym,Zm,metal_in,0.5));
                set(hi,'EdgeColor','none','FaceAlpha',0.65,'FaceColor',[0.90 0.20 0.20]);
                V0 = hi.Vertices; Vras = toRAS(V0); hi.Vertices = Vras(:,1:3);
                legend_handles(end+1) = hi;
                legend_labels{end+1} = 'Metal lead in';
            end
            if any(metal_out(:))
                ho = patch(isosurface(Xm,Ym,Zm,metal_out,0.5));
                set(ho,'EdgeColor','none','FaceAlpha',0.60,'FaceColor',[0.10 0.45 0.85]);
                V0 = ho.Vertices; Vras = toRAS(V0); ho.Vertices = Vras(:,1:3);
                legend_handles(end+1) = ho;
                legend_labels{end+1} = 'Metal lead out';
            end
            axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
            xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
            title('Skull (tan), Metal inside (red), Metal outside (blue)'); view(135,20);
            legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 10);
            filename = [opt.subj_id, '_parts_of_lead_inside_or_ouside_skull_lead_', opt.traj_id];
            full_path = fullfile(opt.report_dir, [filename, '.pdf']);
            saveas(gcf, full_path, 'pdf');
        end
    
        if opt.VisualizeFinal
            legend_handles = [];
            legend_labels = {};
            figure('Color','w'); hold on;
            if any(BW_skull(:))
                hs = patch(isosurface(Xm,Ym,Zm,BW_skull,0.5));
                set(hs,'EdgeColor','none','FaceColor',[0.85 0.80 0.70],'FaceAlpha',0.30);
                V0 = hs.Vertices; Vras = toRAS(V0); hs.Vertices = Vras(:,1:3);
                legend_handles(end+1) = hs;
                legend_labels{end+1} = 'Skull';
            end
            if any(metal_in(:))
                hm = patch(isosurface(Xm,Ym,Zm,metal_in,0.5));
                set(hm,'EdgeColor','none','FaceColor',[0.90 0.20 0.20],'FaceAlpha',0.55);
                V0 = hm.Vertices; Vras = toRAS(V0); hm.Vertices = Vras(:,1:3);
                legend_handles(end+1) = hm;
                legend_labels{end+1} = 'lead inside';
                
            end
            % fitted segment
            h_line = plot3([proximal_entry_ras(1) distal_ras(1)], ...
                  [proximal_entry_ras(2) distal_ras(2)], ...
                  [proximal_entry_ras(3) distal_ras(3)], 'k-', 'LineWidth',2);
            legend_handles(end+1) = h_line;
            legend_labels{end+1} = 'Fitted Trajectory';
            % arrow & spheres
            qScale = 5;
            h_arrow = quiver3(proximal_entry_ras(1), proximal_entry_ras(2), proximal_entry_ras(3), e1(1),e1(2),e1(3), qScale, 'k', 'LineWidth',1.5);
            legend_handles(end+1) = h_arrow;
            legend_labels{end+1} = 'Direction Vector';
            [ssx, ssy, ssz] = sphere(24); r_sph = 2;
            h_prox = surf(ssx*r_sph + proximal_entry_ras(1), ssy*r_sph + proximal_entry_ras(2), ssz*r_sph + proximal_entry_ras(3), ...
                 'EdgeColor','none','FaceColor',[0.10 0.70 0.10],'FaceAlpha',0.95);
            legend_handles(end+1) = h_prox;
            legend_labels{end+1} = 'Proximal Point';
            h_distal = surf(ssx*r_sph + distal_ras(1), ssy*r_sph + distal_ras(2), ssz*r_sph + distal_ras(3), ...
                 'EdgeColor','none','FaceColor',[0.70 0.10 0.70],'FaceAlpha',0.95);
            legend_handles(end+1) = h_distal;
            legend_labels{end+1} = 'Distal Point';
            axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
            xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
            title('Skull (tan), Metal-in (red), Fitted segment (black), Prox (green), Distal (purple)');
            view(135,20);
            legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 10);
            filename = [opt.subj_id, '_parts_of_lead_' , opt.traj_id, '_inside_or_skull_proximal_distal'];
            full_path = fullfile(opt.report_dir, [filename, '.pdf']);
            saveas(gcf, full_path, 'pdf');
        end
        

    end
    
    % possibility for adding a third pt and use spline intepolation,
    % see if spline interpolation is better than linear fitting (could
    % determine by baysien, like information criterion, is it worthy to add one parameter)

    E = proximal_entry_ras; 
    D = distal_ras;
    
    % Choose Rmax sensibly; 1.5–2.0× your tube radius works well
    Rmax_mm = min(2 * 0.8, rmm);
    tube_auc_line = tube_auc_score(P, E, D, 'Rmax', Rmax_mm, 'TrimPct', 5);
    
    % Also gather an RSS for line (for BIC)
    d2_line = dist2_point_line_segment(P, E, D);
    RSS_line = robust_rss(d2_line, 5);
    N = size(P,1);
    BIC_line = N * log(max(RSS_line/N, eps));   % k_line = 0 (E,D fixed by the fit)
    
    % --- 5b) Fit 3-point spline (quadratic Bézier) and score it
    [Ctrl_bez, RSS_bez] = fit_bezier_control_point(P, E, D, ...
        'SamplesM', 201, 'TrimPct', 5, 'Lambda', 0.0);
    
    tube_auc_bez = tube_auc_score(P, E, D, ...
        'BezierC', Ctrl_bez, 'Rmax', Rmax_mm, 'TrimPct', 5);
    
    % BIC for Bézier: 3 params for control point (Cx,Cy,Cz)
    BIC_bez = N * log(max(RSS_bez/N, eps)) + 3 * log(N);
    DeltaBIC = BIC_line - BIC_bez;   % >0 favors Bézier

    % viz tube fitting
    if opt.VisualizeTube_rmm
        nTheta = 40; nSeg = 40; % mesh resolution
        figure('Color','w'); hold on;
        [Eseg, Dseg] = deal(proximal_entry_ras, distal_ras);
        [Vt, Ft] = tube_mesh_segment(Eseg, Dseg, rmm, nTheta, nSeg);
        if ~isempty(Vt)
            htube = patch('Vertices', Vt, 'Faces', Ft, ...
                'FaceColor', [0.10 0.45 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.18);
        end
        
        % (Optional) visualize which metal-in points fall inside the tube
        d2_line_vis = dist2_point_line_segment(P, Eseg, Dseg);
        inTube = sqrt(d2_line_vis) <= rmm;
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
            [tB, Bcen] = bezier_samples(Eseg, Ctrl_bez, Dseg, 60); % get from earlier helper
            % Draw short tubes along each small segment (appears continuous)
            for k = 1:size(Bcen,1)-1
                [Vb, Fb] = tube_mesh_segment(Bcen(k,:), Bcen(k+1,:), rmm, 24, 1);
                if ~isempty(Vb)
                    patch('Vertices', Vb, 'Faces', Fb, ...
                          'FaceColor', [0.10 0.85 0.45], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
                end
            end
            plot3(Bcen(:,1),Bcen(:,2),Bcen(:,3),'g-','LineWidth',2);
        end

        axis equal vis3d; camlight headlight; lighting gouraud;
        xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
        title(sprintf('Metal-in (red), Line (black), Tube r=%.2f mm (blue), In-tube points (green)', rmm));
        view(135,20);
        filename = [opt.subj_id, '_modeling_lead_' , opt.traj_id, '_by_tube_fiting'];
        full_path = fullfile(opt.report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');
    end

    % diagnostic for AUC calculation
    trimPct = 5;
    R = linspace(0, Rmax_mm, 256);
    
    % Distances (line)
    d_line = sqrt(dist2_point_line_segment(P, E, D));
    d_line = sort(d_line);
    n = numel(d_line);
    k0 = floor(n*trimPct/100)+1; k1 = n - floor(n*trimPct/100);
    d_line = d_line(max(1,k0):max(k0,k1));
    
    % Distances (Bézier)
    d_bez = sqrt(dist2_point_bezier(P, E, Ctrl_bez, D, 201));
    d_bez = sort(d_bez);
    n2 = numel(d_bez);
    k0b = floor(n2*trimPct/100)+1; k1b = n2 - floor(n2*trimPct/100);
    d_bez = d_bez(max(1,k0b):max(k0b,k1b));
    
    % Empirical CDFs
    F_line = arrayfun(@(r) mean(d_line <= r), R);
    F_bez  = arrayfun(@(r) mean(d_bez  <= r), R);
    
    % AUCs (normalized to [0,1]) -- should match tube_auc_* up to quadrature error
    AUC_line_curve = trapz(R, F_line) / Rmax_mm;
    AUC_bez_curve  = trapz(R, F_bez)  / Rmax_mm;
    
    % Plot
    fig_auc = figure('Color','w'); hold on;
    plot(R, F_line, 'k-', 'LineWidth', 2);
    plot(R, F_bez,  'g-', 'LineWidth', 2);
    grid on; box on;
    xlabel('Radius r (mm)'); ylabel('Fraction within r (CDF)');
    title('Tube AUC curves (trimmed) — Line vs Bézier');
    
    [bestBIC, bestIdx] = min([BIC_line, BIC_bez]);
    dBIC_line = BIC_line - bestBIC;
    dBIC_bez  = BIC_bez  - bestBIC;
    winner = ["Line","Bézier"];
    winner = winner(bestIdx);
    
    % Kass & Raftery evidence labels
    % bic_label = @(d) ternary(d < 2, "≈ none", ternary(d < 6, "positive", ...
    %                     ternary(d < 10, "strong", "very strong")));
    
    % (Re)plot legend so winner is bolded
    leg = legend({'Line','Bézier'}, 'Location','southeast');
    set(leg, 'AutoUpdate','off');
    
    % Text annotations
    txty = 0.15;
    text(0.60*Rmax_mm, txty+0.16, sprintf('Winner: %s (lower BIC)', winner), ...
        'FontWeight','bold', 'Color','k');
    text(0.60*Rmax_mm, txty+0.10, sprintf('Line   AUC=%.3f  BIC=%.1f  ΔBIC=%.1f ', ...
        AUC_line_curve, BIC_line, dBIC_line), 'Color','k');
    text(0.60*Rmax_mm, txty+0.04, sprintf('Bézier AUC=%.3f  BIC=%.1f  ΔBIC=%.1f ', ...
        AUC_bez_curve,  BIC_bez,  dBIC_bez),  'Color','g');
    
    % Optional: star the winner curve
    if bestIdx == 1
        plot(NaN,NaN,'k*','MarkerSize',10,'DisplayName','Winner');
    else
        plot(NaN,NaN,'g*','MarkerSize',10,'DisplayName','Winner');
    end
    legend('show');

    use_bezier = (DeltaBIC > 10);  % "strong" evidence threshold; adjust if desired
    
    % Generate centerline samples for record/plotting, sampled pts for the
    % chosen model
    t_samp = linspace(0,1,200).';
    if use_bezier
        [~, centerline_pts] = bezier_samples(E, Ctrl_bez, D, numel(t_samp));
        model_chosen = 'bezier';
    else
        centerline_pts = E + t_samp.*(D - E);
        model_chosen = 'line';
    end
    filename = [opt.subj_id, '_tube_auc_score_changing_radius_for_choosing_line_and_bezier_lead_', opt.traj_id];
    full_path = fullfile(opt.report_dir, [filename, '.pdf']);
    saveas(gcf, full_path, 'pdf');


    % ---- 8) Pack outputs once ----
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
        'boundary_hits',      struct('B1',B1,'B2',B2,'t1',t1,'t2',t2), ...
        'tube_auc_line',     AUC_line_curve, ...
        'tube_auc_bezier',   AUC_bez_curve, ...
        'bic_line',          BIC_line, ...
        'bic_bezier',        BIC_bez, ...
        'delta_bic',         DeltaBIC, ...
        'model',             model_chosen, ...          % 'line' or 'bezier'
        'bezier_control',    Ctrl_bez, ...              % 1x3 (empty if line chosen)
        'centerline_t',      t_samp, ...                % Nx1 in [0,1]
        'centerline_pts',    centerline_pts ...        % Nx3 RAS samples of chosen model
    );
end