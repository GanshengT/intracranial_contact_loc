function [res, stats] = fit_axis_cylinder_robust(P, varargin)
% Robust line (axis) fit to a blob point cloud in RAS.
% Minimizes orthogonal distance with Tukey biweight via IRLS; picks tube radius r.
%
% res: struct with p0 (RAS 1x3), u (unit 1x3), r (mm), inliers (Nx1 logical)

ip = inputParser;
ip.addParameter('MaxIter', 20);
ip.addParameter('Tol', 1e-5);
ip.addParameter('RadiusInit', []);
ip.addParameter('RadiusGrid', linspace(0.5, 2.0, 6)); % mm
ip.parse(varargin{:});
opt = ip.Results;

N = size(P,1);
% --- init via SVD (more stable than eig on cov when elongated)
mu = mean(P,1);
[U,~,Vt] = svd(bsxfun(@minus,P,mu),'econ'); %
u = Vt(:,1)'; u = u./max(norm(u),eps);
p0 = mu;  % for infinite line, centroid is optimal intercept

last_obj = Inf;
for it = 1:opt.MaxIter
    % orthogonal distances to current line
    V = P - p0;          % Nx3
    proj = (V * u.');    % Nx1
    par  = proj .* u;    % Nx3
    R    = V - par;      % residuals (perp)
    d    = sqrt(sum(R.^2,2));  % Nx1

    % robust scale and weights (Tukey biweight)
    sigma = 1.4826 * mad(d,1) + eps;
    c = 4.685 * sigma;
    r = d ./ c;
    w = (abs(r) < 1) .* (1 - r.^2).^2;  % Tukey

    if ~any(w)  % degenerate
        w = ones(N,1);
    end

    % weighted centroid
    Wsum = sum(w) + eps;
    mu_w = (w.'*P) ./ Wsum;

    % weighted covariance (centered)
    Pc   = bsxfun(@minus,P,mu_w);
    Cw   = (Pc' * (bsxfun(@times, Pc, w))) / Wsum;

    % principal axis from weighted cov
    [Vax,Dax] = eig((Cw + Cw')/2);
    [~,ix] = max(diag(Dax));
    u_new = Vax(:,ix)'; u_new = u_new./max(norm(u_new),eps);

    % project weighted centroid onto new axis for p0
    p0_new = mu_w + ((mu_w - mu_w) * u_new.') * u_new; % == mu_w (kept explicit)

    % objective = weighted sum of squared perp distances
    obj = sum(w .* d.^2);

    % update
    if abs(last_obj - obj) < opt.Tol*max(1,last_obj), break; end
    last_obj = obj; u = u_new; p0 = p0_new;
end

% get tube radius that best explains the blob
V  = P - p0; proj = (V*u.'); par = proj.*u; d = sqrt(sum((V-par).^2,2));
bestScore = -Inf; r_best = [];
for r_try = opt.RadiusGrid
    in_try  = d <= r_try;
    if ~any(in_try), continue; end
    t_try   = proj(in_try);
    span    = prctile(t_try,97.5) - prctile(t_try,2.5);
    cover   = mean(in_try);
    score   = cover * span;           % simple tradeoff: coverage × axial span
    if score > bestScore
        bestScore = score; r_best = r_try; in_best = in_try; span_best = span;
    end
end

% if no radius grid hits, set radius to median distance
if isempty(r_best)
    r_best = median(d) + eps; in_best = d <= r_best; span_best = range(proj(in_best));
end

res = struct('p0',p0,'u',u,'r',r_best,'inliers',in_best);
stats = struct('iters',it,'obj',last_obj,'span',span_best,'coverage',mean(in_best));
end