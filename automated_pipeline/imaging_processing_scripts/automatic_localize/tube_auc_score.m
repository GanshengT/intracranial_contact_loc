function s = tube_auc_score(P, E, D, varargin)
% TUBE_AUC_SCORE  Radius-agnostic tube fit score in [0,1].
% s = 1 - mean(min(d, Rmax))/Rmax, where d are point-to-centerline distances.
% Inputs:
%   P : Nx3 RAS points (metal-in)
%   E,D : 1x3 endpoints (proximal, distal)
% Options:
%   'BezierC'   : []  (1x3 control point -> use quadratic Bézier if given)
%   'SamplesM'  : 201 (Bézier sampling density)
%   'Rmax'      : 3.0 (mm)  % pick > expected tube radius; insensitive within reason
%   'TrimPct'   : 5    % trim far outliers before averaging (robust)
%   'WeightFn'  : []   % optional @(r) weight over radius; default = uniform
%
% Output:
%   s : scalar in [0,1]; higher is better fit.

p = inputParser;
addParameter(p,'BezierC',[]);
addParameter(p,'SamplesM',201);
addParameter(p,'Rmax',3.0);
addParameter(p,'TrimPct',5);
addParameter(p,'WeightFn',[]);
parse(p,varargin{:});
C = p.Results.BezierC; M = p.Results.SamplesM;
Rmax = p.Results.Rmax; trim = p.Results.TrimPct;
wfun = p.Results.WeightFn;

% distances to chosen centerline
if isempty(C)
    d2 = dist2_point_line_segment(P, E, D);
else
    d2 = dist2_point_bezier(P, E, C, D, M);
end
d = sqrt(d2);

% robust trimming (drop largest trim% and smallest trim% if desired)
d = sort(d);
n = numel(d);
k0 = floor(n*trim/100)+1;
k1 = n - floor(n*trim/100);
d = d(max(1,k0):max(k0,k1));

if isempty(wfun)
    % Uniform weighting over radius (closed form)
    s = 1 - mean(min(d, Rmax)) / Rmax;
else
    % Optional custom weighting via numerical quadrature over r
    % (e.g., w(r)=exp(-r/σ) to emphasize near-axis fit)
    R = linspace(0, Rmax, 256);
    F = arrayfun(@(r) mean(d <= r), R);
    W = arrayfun(wfun, R);
    s = trapz(R, F.*W) / trapz(R, W);      % normalized weighted AUC
end
end

% ---- helpers reused ----
function d2 = dist2_point_line_segment(P, E, D)
U = D - E; L2 = sum(U.^2);
t = ((P - E) * U.') / max(L2, eps);
t = min(max(t,0),1);
Proj = E + t.*U;
d2 = sum((P - Proj).^2, 2);
end

function [t,B] = bezier_samples(E, C, D, M)
t = linspace(0,1,M).';
B = (1-t).^2 .* E + 2*(1-t).*t .* C + t.^2 .* D;
end

function d2 = dist2_point_bezier(P, E, C, D, M)
[~,B] = bezier_samples(E,C,D,M);
d2 = min(pdist2(P, B, 'euclidean').^2, [], 2);
end