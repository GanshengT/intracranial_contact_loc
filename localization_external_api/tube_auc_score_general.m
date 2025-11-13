function auc = tube_auc_score_general(P, model, args, Rmax, TrimPct)
% Unified Tube AUC scorer for line, Bézier, or polyline models
if nargin<5, TrimPct = 5; end

% Compute distances
switch lower(model)
    case 'line'
        d = sqrt(dist2_point_line_segment(P, args.E, args.D));
    case 'bezier'
        M = 201; 
        if isfield(args,'SamplesM'), M = args.SamplesM; end
        d = sqrt(dist2_point_bezier(P, args.E, args.C, args.D, M));
    case 'polyline'
        d = sqrt(dist2_point_polyline(P, args.C));
    otherwise
        error('Unknown model: %s', model);
end

% Trim extremes and compute empirical CDF
d = sort(d(:));
n = numel(d);
k0 = floor(n*TrimPct/100)+1; 
k1 = n - floor(n*TrimPct/100);
if k0 > k1, k0 = 1; k1 = n; end
d = d(k0:k1);

R = linspace(0, Rmax, 256);
F = arrayfun(@(r) mean(d <= r), R);

% Normalized AUC
auc = trapz(R, F) / max(Rmax, eps);
end