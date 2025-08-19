function [t0, t1, r_mm, stats] = enclosing_tube_params(G, Xras, ctr, dstar, varargin)
% ENCLOSING_TUBE_PARAMS  Tight tube around the current blob along dstar.
% Returns longitudinal span [t0,t1] in mm (relative to ctr) and radius r_mm.
%
% Name-Value (optional):
%   'TGuardMM' (0.5) : pad added to each end of [t0,t1]
%   'RGuardMM' (0.2) : pad added to radius
%   'RobustPct' (100): use max radial (100) or a percentile (e.g., 98) to ignore outliers

    ip = inputParser; ip.CaseSensitive=false;
    addParameter(ip,'TGuardMM',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'RGuardMM',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'RobustPct',100,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=100);
    parse(ip,varargin{:});
    TGuard = ip.Results.TGuardMM;
    RGuard = ip.Results.RGuardMM;
    q      = ip.Results.RobustPct;

    if ~any(G(:))
        t0 = -1; t1 = 1; r_mm = 1; stats = struct('n',0,'tmin',t0,'tmax',t1,'rmax',r_mm); 
        return;
    end

    % collect blob points in RAS
    [ii,jj,kk] = ind2sub(size(G), find(G));
    P = [ squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*1))), ...
          squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*2))), ...
          squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*3))) ];

    d = dstar(:).'/max(norm(dstar),eps);
    % shift to axis origin at ctr
    X = P - ctr;
    % longitudinal coordinates
    t = X*d(:);                    % Nx1
    % radial distances to axis
    cx = cross(X, repmat(d,numel(t),1), 2);
    r  = sqrt(sum(cx.^2,2));       % Nx1

    t0 = min(t) - TGuard;
    t1 = max(t) + TGuard;

    if q >= 100
        r_mm = max(r) + RGuard;
    else
        r_mm = prctile(r, q) + RGuard;
    end

    stats = struct('n',size(P,1),'tmin',min(t),'tmax',max(t), ...
                   'rmax',max(r),'rq',prctile(r,q));
end

