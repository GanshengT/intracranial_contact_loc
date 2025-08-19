function [dstar, bestFill, t0_best, t1_best, r_best] = fit_axis_by_fill( ...
        U_in, ctr, TPadMM, G, Xras, voxel_mm, varargin)
% FIT_AXIS_BY_FILL  (volume-fill; integer-degree ring search)
% Searches a cone around U_in by sweeping rings at 1°..AngleDeg (inclusive).
% For each direction, builds a tight tube around the blob and maximizes:
%     fill = V_G / V_tube.
%
% Inputs:
%   U_in      : 1x3 planned axis (RAS)
%   ctr       : 1x3 current blob centroid (RAS, mm)
%   TPadMM    : pad added to t-extent (mm)
%   G         : logical blob mask (MxNxP)
%   Xras      : MxNxP x 3 voxel centers (mm)
%   voxel_mm  : [dI dJ dK] voxel sizes (mm) in array order
%
% Name-Value:
%   'AngleDeg'   (5)    : half-cone (deg). Rings at 1..AngleDeg
%   'RingK'      (12)   : azimuthal samples per ring
%   'IncludeCenter' (true) : evaluate U_in itself as well
%   'RGuardMM'   (0.2)  : radial guard (mm)
%   'TGuardMM'   (0.5)  : extra end guard (mm) in addition to TPadMM
%   'RobustPct'  (100)  : radius from percentile of r (e.g., 98)
%
% Outputs:
%   dstar, bestFill, t0_best, t1_best, r_best

    ip = inputParser; ip.CaseSensitive=false;
    addParameter(ip,'AngleDeg',5,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=45);
    addParameter(ip,'RingK',12,@(x)isnumeric(x)&&isscalar(x)&&x>=4&&x<=128);
    addParameter(ip,'IncludeCenter',true,@(x)islogical(x)||ismember(x,[0 1]));
    addParameter(ip,'RGuardMM',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'TGuardMM',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'RobustPct',100,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=100);
    parse(ip,varargin{:});
    AngleDeg   = ip.Results.AngleDeg;
    RingK      = ip.Results.RingK;
    IncCenter  = logical(ip.Results.IncludeCenter);
    RGuardMM   = ip.Results.RGuardMM;
    TGuard2    = ip.Results.TGuardMM;
    qPct       = ip.Results.RobustPct;

    EPS_NUM = 1e-9;

    % quick exit if G empty
    if ~any(G(:))
        U0 = U_in(:).'/max(norm(U_in),EPS_NUM);
        dstar=U0; bestFill=0; t0_best=-1; t1_best=1; r_best=1;
        return;
    end

    % blob points in RAS and relative to ctr
    [ii,jj,kk] = ind2sub(size(G), find(G));
    PG = [ squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*1))), ...
           squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*2))), ...
           squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*3))) ];
    X = PG - ctr;                         % Nx3
    NG = size(PG,1);
    Vvoxel = prod(voxel_mm);              % mm^3
    VG = NG * Vvoxel;                     % mm^3

    % planned axis (unit) and local frame
    U0 = U_in(:).'/max(norm(U_in),EPS_NUM);
    tmp = [1 0 0]; if abs(dot(tmp,U0))>0.9, tmp=[0 1 0]; end
    v1 = tmp - dot(tmp,U0)*U0; v1 = v1 / norm(v1);
    v2 = cross(U0, v1);

    % build candidate directions: center + rings at 1..AngleDeg (deg)
    dirs = zeros(0,3);
    if IncCenter
        dirs(end+1,:) = U0; 
    end
    if AngleDeg > 0
        az = linspace(0, 2*pi, RingK+1); az(end) = [];
        for deg = 1:AngleDeg
            aa = deg2rad(deg);
            ca = cos(aa); sa = sin(aa);
            for th = az
                d = ca*U0 + sa*(cos(th)*v1 + sin(th)*v2);
                d = d ./ max(norm(d), EPS_NUM);
                dirs(end+1,:) = d;
            end
        end
    end

    % evaluate volume-fill for each direction
    bestFill = -inf; dstar = U0; t0_best = 0; t1_best = 0; r_best = 1;

    for k = 1:size(dirs,1)
        d = dirs(k,:);

        % t-extent (with padding)
        t = X * d(:);                                 % NG x 1
        t0 = min(t) - (TPadMM + TGuard2);
        t1 = max(t) + (TPadMM + TGuard2);
        L  = max(t1 - t0, EPS_NUM);

        % radial distances to axis
        cx = cross(X, repmat(d,NG,1), 2);
        r  = sqrt(sum(cx.^2, 2));
        if qPct < 100
            r_use = prctile(r, qPct);
        else
            r_use = max(r);
        end
        r_mm = max(r_use + RGuardMM, EPS_NUM);

        % tube volume and fill
        Vtube = pi * (r_mm^2) * L;
        f = VG / Vtube;

        if f > bestFill
            bestFill = f;
            dstar    = d;
            t0_best  = t0;
            t1_best  = t1;
            r_best   = r_mm;
        end
    end
end