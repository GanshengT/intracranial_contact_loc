function [C_hat, RSS_bez] = fit_bezier_control_point(P, E, D, varargin)
% Optimize quadratic Bézier control point C (E and D fixed).
    p = inputParser;
    addParameter(p,'SamplesM',201);
    addParameter(p,'TrimPct',5);
    addParameter(p,'Lambda',0.0);  % curvature regularization
    parse(p,varargin{:});
    M = p.Results.SamplesM; trim = p.Results.TrimPct; lambda = p.Results.Lambda;
    
    U = D - E; Un = U / max(norm(U), eps);
    tmp = [1 0 0]; if abs(dot(tmp,Un))>0.9, tmp=[0 1 0]; end
    N1 = cross(Un, tmp); N1 = N1 / max(norm(N1), eps);
    mid = 0.5*(E + D);
    C0  = mid + 2*N1;  % ~2 mm off-axis
    
    obj = @(C) robust_rss(dist2_point_bezier(P, E, C, D, M), trim) ...
               + lambda * sum((E - 2*C + D).^2);
    
    C_hat = fminsearch(obj, C0, optimset('Display','off','TolX',1e-4,'TolFun',1e-4));
    RSS_bez = robust_rss(dist2_point_bezier(P, E, C_hat, D, M), trim);
end