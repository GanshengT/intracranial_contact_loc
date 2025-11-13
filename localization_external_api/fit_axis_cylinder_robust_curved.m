function [res, stats] = fit_axis_cylinder_robust_curved(P, varargin)
% Robust curved tube fit to a voxel cloud (RAS).
% Principal-curve style via IRLS with Tukey biweight on orthogonal distance.
% Adds arc-length reparameterization and turning-aware ("knee-boosting") weights.
%
% Inputs
%   P: Nx3 voxel coordinates (RAS)
%
% Options (name/value)
%   'MaxIter'          (20)        IRLS outer iterations
%   'Tol'              (1e-4)      Relative change tolerance on robust objective
%   'SampleM'          (1000)      Dense samples used for projection speedup
%   'SplineSmoothing'  (1e-3)      csaps smoothing parameter (0..1)
%   'KnotCount'        ([])        reserved, not used (csaps sets knots)
%   'RadiusMode'       ('percentile')  'grid'|'percentile'
%   'RadiusGrid'       (linspace(0.3,2.0,9))  used if RadiusMode='grid'
%   'RadiusPercentile' (90)        percentile for radius estimate
%   'VarRadius'        (true)      estimate r(s) locally (moving percentile)
%   'VarRadiusWin'     (0.05)      window width over s in [0,1] for local r
%   'CosineFilter'     (false)     prefilter voxels by directional coherence
%   'CosinePercent'    (30)        % farthest neighbors to compute coherence
%   'CosineThresh'     (0.97)      keep voxels with coherence >= this
%   'Report'           ([])        struct with subj_id, report_dir, traj_id
%   'SaveFigures'      (true)      save non-interactive figures if Report provided
%   'KneeBoost'        (true)      enable turning-aware weight boost
%
% Outputs
%   res:
%     .ppx, .ppy, .ppz    ppform splines x(s),y(s),z(s) with s in [0,1]
%     .C                  Mx3 centerline samples at s_dense
%     .T                  Mx3 unit tangents at s_dense
%     .kappa              Mx1 curvature at s_dense (1/mm)
%     .t_proj             Nx1 param in [0,1] of closest point on curve
%     .dist_perp          Nx1 orthogonal distances (mm)
%     .r                  scalar radius (if VarRadius=false)
%     .r_of_t             Mx1 local radius profile (if VarRadius=true)
%     .inliers            Nx1 logical (dist<= local/global radius)
%     .knees_idx          indices of detected knee peaks along s_dense
%   stats:
%     .iters             IRLS iters
%     .obj               final robust objective
%     .sigma             robust scale of distances
%     .converged         logical
%
% G.T. g.tan@wustl.edu

% ----- options -----
ip = inputParser;
ip.addParameter('MaxIter', 20);
ip.addParameter('Tol', 1e-4);
ip.addParameter('SampleM', 1000);
ip.addParameter('SplineSmoothing', 1e-3);
ip.addParameter('KnotCount', []); % (kept for API harmony)
ip.addParameter('RadiusMode','percentile'); % 'grid' or 'percentile'
ip.addParameter('RadiusGrid', linspace(0.3,2.0,9));
ip.addParameter('RadiusPercentile', 90);
ip.addParameter('VarRadius', true);
ip.addParameter('VarRadiusWin', 0.05);
ip.addParameter('CosineFilter', false);
ip.addParameter('CosinePercent', 30);
ip.addParameter('CosineThresh', 0.97);
ip.addParameter('Report', []);
ip.addParameter('SaveFigures', true);
ip.addParameter('KneeBoost', true);
ip.parse(varargin{:});
opt = ip.Results;

% ----- input -----
P = double(P);
N = size(P,1);
if N < 20, error('Not enough points to fit a curve.'); end

% ----- optional prefilter: directional coherence -----
if opt.CosineFilter
    cos_stats = local_directional_coherence(P, opt.CosinePercent);
    keep_mask = cos_stats >= opt.CosineThresh;
    if any(~keep_mask)
        P = P(keep_mask,:);
        N = size(P,1);
    end
    if N < 20, error('Too few points after coherence filtering.'); end
end

% ----- initialization: robust straight line via SVD -----
mu = mean(P,1);
[~,~,Vt] = svd(bsxfun(@minus,P,mu),'econ');
u0 = Vt(:,1)';  u0 = u0 ./ max(norm(u0),eps);
t0 = (P - mu) * u0.';                 % scalar coordinate along initial line
t0 = normalize_to_unit_interval(t0);  % map to [0,1]

% initial csaps with equal weights (param = t0)
w = ones(N,1);
ppx = csaps(t0.', P(:,1).', opt.SplineSmoothing, [], w.');
ppy = csaps(t0.', P(:,2).', opt.SplineSmoothing, [], w.');
ppz = csaps(t0.', P(:,3).', opt.SplineSmoothing, [], w.');

% ----- one-time "before IRLS" visualization (non-interactive) -----
if save_ok(opt)
    try
        Md = max(200, round(opt.SampleM));
        td = linspace(0,1,Md);
        Cx = ppval(ppx, td(:)); Cy = ppval(ppy, td(:)); Cz = ppval(ppz, td(:));
        Cinit = [Cx(:), Cy(:), Cz(:)];
        fig = figure('Color','w'); hold on; axis equal vis3d;
        scatter3(P(:,1),P(:,2),P(:,3),4,[0.7 0.7 0.7],'filled','MarkerFaceAlpha',0.2);
        plot3(Cinit(:,1),Cinit(:,2),Cinit(:,3),'k-','LineWidth',1.5);
        title('Initialization: points (gray) and initial spline (black)');
        xlabel X; ylabel Y; zlabel Z; view(135,20); grid on;
        fn = sprintf('%s_curve_init_%s', opt.Report.subj_id, opt.Report.traj_id);
        saveas(fig, fullfile(opt.Report.report_dir, [fn '.pdf']), 'pdf');
        close(fig);
    catch, end
end

% ----- IRLS principal-curve loop -----
last_obj = Inf;
converged = false;

for it = 1:opt.MaxIter
    % (1) project points to current curve param in [0,1]
    t_i = project_points_to_curve(P, ppx, ppy, ppz, opt.SampleM); % N×1

    % (2) orthogonal distances & Tukey weights
    Cx = ppval(ppx, t_i(:)); Cy = ppval(ppy, t_i(:)); Cz = ppval(ppz, t_i(:));
    C_i = [Cx(:), Cy(:), Cz(:)];               % Nx3
    R_i = P - C_i;     
    d_i = sqrt(sum(R_i.^2,2));                 % Nx1

    sigma = 1.4826 * mad(d_i,1) + eps;
    cscale = 4.685 * sigma;
    rr = d_i ./ cscale;
    w = (abs(rr) < 1) .* (1 - rr.^2).^2;       % Tukey biweight
    if ~any(w), w = ones(N,1); end

    % (2.5) turning-aware knee boost (optional)
    if opt.KneeBoost
        [boost, ~] = turning_weight_boost(t_i, C_i);
        % clamp for stability
        boost = max(1, min(3, boost));
        w = w .* boost;
    end

    % robust objective (for convergence)
    obj = sum(w .* (d_i.^2));

    % (3) refit splines x(t),y(t),z(t) with weights at (t_i, P)
    ppx = csaps(t_i.', P(:,1).', opt.SplineSmoothing, [], w.');
    ppy = csaps(t_i.', P(:,2).', opt.SplineSmoothing, [], w.');
    ppz = csaps(t_i.', P(:,3).', opt.SplineSmoothing, [], w.');

    % (3.5) arc-length reparameterization to s in [0,1]
    Mdense = max(400, round(1.5*opt.SampleM));
    t_tmp  = linspace(0,1,Mdense);
    CxT = ppval(ppx, t_tmp(:)); CyT = ppval(ppy, t_tmp(:)); CzT = ppval(ppz, t_tmp(:));
    Ctmp = [CxT(:), CyT(:), CzT(:)];
    seg = sqrt(sum(diff(Ctmp,1,1).^2, 2));
    S   = [0; cumsum(seg)];
    if S(end) > 0, S = S ./ S(end); else, S = linspace(0,1,numel(t_tmp)).'; end

    % map t_i -> s_i, then refit against s_i with same weights
    s_i = interp1(t_tmp, S, t_i, 'linear', 'extrap');
    s_i = max(0, min(1, s_i));
    ppx = csaps(s_i.', P(:,1).', opt.SplineSmoothing, [], w.');
    ppy = csaps(s_i.', P(:,2).', opt.SplineSmoothing, [], w.');
    ppz = csaps(s_i.', P(:,3).', opt.SplineSmoothing, [], w.');

    % (4) convergence
    if abs(last_obj - obj) < opt.Tol*max(1,last_obj)
        converged = true; break;
    end
    last_obj = obj;
end

% ----- final projection & distances (now param is s in [0,1]) -----
t_i = project_points_to_curve(P, ppx, ppy, ppz, opt.SampleM);
Cx = ppval(ppx, t_i(:)); Cy = ppval(ppy, t_i(:)); Cz = ppval(ppz, t_i(:));
C_i = [Cx(:), Cy(:), Cz(:)];
R_i = P - C_i;
d_i = sqrt(sum(R_i.^2,2));

% ----- dense centerline samples & differential geometry -----
M = max(200, min(2000, round(opt.SampleM)));
s_dense = linspace(0,1,M);
Cx = ppval(ppx, s_dense(:)); Cy = ppval(ppy, s_dense(:)); Cz = ppval(ppz, s_dense(:));
C  = [Cx(:), Cy(:), Cz(:)];

dppx = fnder(ppx,1); dppy = fnder(ppy,1); dppz = fnder(ppz,1);
ddppx = fnder(ppx,2); ddppy = fnder(ppy,2); ddppz = fnder(ppz,2);

d1 = [ppval(dppx,s_dense(:)), ppval(dppy,s_dense(:)), ppval(dppz,s_dense(:))];
d2 = [ppval(ddppx,s_dense(:)), ppval(ddppy,s_dense(:)), ppval(ddppz,s_dense(:))];

T = normalize_rows(d1);
cross12 = cross(d1, d2, 2);
kappa = vecnorm(cross12,2,2) ./ max(vecnorm(d1,2,2).^3, eps);

% ----- knee detection (robust) -----
thr_pct = prctile(kappa, 95);
thr_abs = 0.01 * mean(kappa + eps);
thr = max(thr_pct, thr_abs);

is_peak = false(numel(kappa),1);
for i = 2:numel(kappa)-1
    if kappa(i) >= kappa(i-1) && kappa(i) >= kappa(i+1) && kappa(i) >= thr
        is_peak(i) = true;
    end
end
knees_idx = find(is_peak);

% ----- post-fit visualizations (non-interactive; saved if Report provided) -----
if save_ok(opt)
    try
        % (A) points + final curve + knees
        fig1 = figure('Color','w'); hold on; axis equal vis3d;
        scatter3(P(:,1),P(:,2),P(:,3),4,[0.75 0.75 0.75],'filled','MarkerFaceAlpha',0.25);
        plot3(C(:,1),C(:,2),C(:,3),'k-','LineWidth',1.8);
        if ~isempty(knees_idx)
            scatter3(C(knees_idx,1),C(knees_idx,2),C(knees_idx,3),40,'m','filled');
        end
        title('Final curved fit with knee markers');
        xlabel X; ylabel Y; zlabel Z; grid on; view(135,20);
        fn1 = sprintf('%s_curve_final_%s', opt.Report.subj_id, opt.Report.traj_id);
        saveas(fig1, fullfile(opt.Report.report_dir, [fn1 '.pdf']), 'pdf');
        close(fig1);

        % (B) curvature profile & histogram
        fig2 = figure('Color','w','Position',[100 100 900 360]);
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

        nexttile; hold on;
        plot(s_dense, kappa, 'k-', 'LineWidth', 1.5);
        yline(thr, '--r', 'Thresh'); 
        if ~isempty(knees_idx)
            scatter(s_dense(knees_idx), kappa(knees_idx), 30, 'm', 'filled');
        end
        xlabel('normalized arc length s'); ylabel('\kappa (1/mm)');
        title('\kappa(s) with knee peaks'); grid on;

        nexttile; hold on;
        histogram(kappa, 40, 'Normalization','pdf');
        xline(thr, '--r', 'Thresh');
        xlabel('\kappa (1/mm)'); ylabel('pdf');
        title('Curvature distribution'); grid on;

        fn2 = sprintf('%s_curvature_profile_%s', opt.Report.subj_id, opt.Report.traj_id);
        saveas(fig2, fullfile(opt.Report.report_dir, [fn2 '.pdf']), 'pdf');
        close(fig2);
    catch
        % best-effort plotting
    end
end

% ----- radius estimation -----
switch lower(opt.RadiusMode)
    case 'grid'
        bestScore = -Inf; r_const = [];
        for r_try = opt.RadiusGrid
            in_try = d_i <= r_try;
            if ~any(in_try), continue; end
            span = prctile(t_i(in_try),97.5) - prctile(t_i(in_try),2.5);
            cover = mean(in_try);
            score = cover * span;
            if score > bestScore
                bestScore = score; r_const = r_try; in_best = in_try;
            end
        end
        if isempty(r_const), r_const = median(d_i); in_best = d_i <= r_const; end
        r_of_t = [];
    otherwise % 'percentile'
        r_const = prctile(d_i, opt.RadiusPercentile);
        in_best = d_i <= r_const;
        r_of_t = [];
end

% variable radius along centerline (local percentile)
if opt.VarRadius
    win = max(1/M, min(0.25, opt.VarRadiusWin));  % window in s-space
    r_of_t = zeros(M,1);
    for m = 1:M
        s0 = s_dense(m);
        mask = abs(t_i - s0) <= win;
        if ~any(mask)
            r_of_t(m) = r_const;
        else
            r_of_t(m) = prctile(d_i(mask), opt.RadiusPercentile);
        end
    end
    [~, idxNear] = min(abs(t_i - s_dense(:).'), [], 2);
    r_local = r_of_t(idxNear);
    in_best = d_i <= r_local;
end

% ----- outputs -----
res = struct();
res.ppx = ppx; res.ppy = ppy; res.ppz = ppz;     % centerline splines x(s),y(s),z(s)
res.C = C; res.T = T; res.kappa = kappa;         % samples / tangents / curvature
res.t_proj = t_i(:); res.dist_perp = d_i(:);
if opt.VarRadius, res.r = [];  res.r_of_t = r_of_t(:);
else,             res.r = r_const; res.r_of_t = [];
end
res.inliers = in_best(:);
res.knees_idx = knees_idx;

stats = struct('iters', it, 'obj', last_obj, 'sigma', sigma, 'converged', converged);

end % ===== main =====


% ===================== helpers =====================

function ok = save_ok(opt)
ok = ~isempty(opt.Report) && isstruct(opt.Report) ...
     && all(isfield(opt.Report, {'subj_id','report_dir','traj_id'})) ...
     && isfield(opt, 'SaveFigures') && logical(opt.SaveFigures);
end

function t_i = project_points_to_curve(P, ppx, ppy, ppz, M)
% Coarse nearest among M dense samples, then 1D local refine (fminbnd)
N = size(P,1);
t_dense = linspace(0,1,M);
Cx = ppval(ppx, t_dense(:)); Cy = ppval(ppy, t_dense(:)); Cz = ppval(ppz, t_dense(:));
C = [Cx(:), Cy(:), Cz(:)];
[~, idx] = pdist2(C, P, 'euclidean', 'Smallest', 1);
t0 = t_dense(idx).';
t_i = zeros(N,1);
for n = 1:N
    tn = t0(n);
    a = max(0, tn - 0.05);
    b = min(1, tn + 0.05);
    f = @(t) sum((P(n,:) - [ppval(ppx,t), ppval(ppy,t), ppval(ppz,t)]).^2, 2);
    t_i(n) = fminbnd(f, a, b);
end
t_i = max(0, min(1, t_i));
end

function t = normalize_to_unit_interval(s)
smin = min(s); smax = max(s);
if smax > smin, t = (s - smin) / (smax - smin);
else,           t = zeros(size(s));
end
end

function U = normalize_rows(V)
n = sqrt(sum(V.^2,2)) + eps;
U = V ./ n;
end

function [boost, ang_s] = turning_weight_boost(t_i, C_i)
% Turning-aware weight multiplier:
% - sort by current param t_i
% - compute finite-diff tangents and angle-of-turn proxy
% - smooth and map back to original ordering
[t_sort, ord] = sort(t_i);
C_sorted = C_i(ord,:);
dC = gradient(C_sorted, max(1e-6, median(diff(t_sort))+eps));
dC = normalize_rows(dC);
ang = zeros(size(dC,1),1);
if size(dC,1) >= 3
    % compare vectors 2 steps apart (more stable than 1-step)
    v1 = dC(1:end-2,:); v2 = dC(3:end,:);
    cosang = sum(v1.*v2,2);
    cosang = max(-1,min(1,cosang));
    ang(2:end-1) = acos(cosang);
end
win = max(5, round(0.01*length(ang)));
ang_s = movmean(ang, win, 'Endpoints','shrink');

ang_full = zeros(size(t_i));
ang_full(ord) = ang_s;

s_ang = median(ang_full) + mad(ang_full,1);
boost = 1 + 2.0 * (ang_full ./ max(s_ang,eps)); % 1x..3x (will clamp upstream)
end

function s = local_directional_coherence(P, percent)
% crude long-range directional coherence: avg pairwise cosine among far neighbors
N = size(P,1);
K = max(10, round(percent/100 * N));
[~, idx] = pdist2(P, P, 'euclidean', 'Smallest', K+1);
idx = idx(2:end,:);
s = zeros(N,1);
for i = 1:N
    V = bsxfun(@minus, P(idx(:,i),:), P(i,:));    % Kx3
    V = V ./ (vecnorm(V,2,2)+eps);
    C = V*V.';                                   % KxK cosine matrix
    up = triu(true(K),1);
    s(i) = mean(C(up));
end
end