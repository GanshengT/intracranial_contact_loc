function distal_ras = distal_from_pointcloud_robust(P, proximal_entry_ras, e1, varargin)
% Robustly pick distal along +e1 using the metal point cloud itself.
% P  : Nx3 metal points (RAS)
% e1 : 1x3 unit axis direction (RAS), oriented proximal->distal
% proximal_entry_ras : RAS of proximal boundary (origin for t)

ip = inputParser;
ip.addParameter('BinMM', 0.25);      % histogram bin width (mm)
ip.addParameter('SmoothWin', 9);     % odd window length for 1D smoothing
ip.addParameter('TailPct', 98.5);    % tail percentile for distal (97.5–99 works well)
ip.addParameter('ClipTMax', inf);    % optional brain boundary t2 to clip distal
ip.addParameter('RadQuantile', 0.90);% define near-axis radius by this quantile
ip.parse(varargin{:});
opt = ip.Results;

% axial coords from proximal, radial distances to axis
V   = P - proximal_entry_ras;         % Nx3
t   = V * e1.';                       % Nx1 (mm), signed along +e1
par = (t .* e1);                      % Nx3 parallel component
rad = vecnorm(V - par, 2, 2);         % Nx1 radial distance (mm)

% robust near-axis radius
r_est = quantile(rad, opt.RadQuantile);
in    = rad <= r_est;                 % keep points near the axis
if ~any(in)
    % fallback: use all points
    in = true(size(t));
end
t_in = t(in);
t_in = t_in(t_in >= 0);               % distal is along +e1 from proximal

if numel(t_in) < 5
    % degenerate fallback: farthest point along +e1
    t_star = max(0, max(t));
else
    % histogram + smoothing
    if isempty(t_in), t_in = 0; end
    edges = (min(t_in)-opt.BinMM) : opt.BinMM : (max(t_in)+opt.BinMM);
    if numel(edges) < 3, edges = linspace(0, max(t_in)+opt.BinMM, 8); end
    cnt = histcounts(t_in, edges);
    ctr = movmean(edges,2,'Endpoints','discard'); % bin centers

    % smooth with boxcar (robust)
    w = opt.SmoothWin;
    if mod(w,2)==0, w = w+1; end
    cnt_s = movmean(cnt, w, 'Endpoints','shrink');

    % tail threshold from percentile on t_in
    t_tail = prctile(t_in, opt.TailPct);

    % find last index where smoothed density is above a small fraction of its peak
    thr = 0.08 * max(cnt_s);   % 8% of peak; tweak 5–15% if needed
    idx = find(ctr >= 0 & cnt_s >= thr);
    if isempty(idx)
        t_star = t_tail;
    else
        t_star = max(t_tail, ctr(idx(end)));
    end
end

%  clip to boundary, optional
t_star = min(t_star, opt.ClipTMax);

distal_ras = proximal_entry_ras + t_star * e1;
end