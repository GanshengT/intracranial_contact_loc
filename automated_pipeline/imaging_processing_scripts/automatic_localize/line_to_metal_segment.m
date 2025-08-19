function [enter_ras, exit_ras] = line_to_metal_segment(BW_metal, Avox2ras, p_start, u, t_range, step_mm)
% Along p(t)=p_start + t*u, find the widest contiguous run inside BW_metal (logical, CT grid).
% Returns RAS points at the run's enter and exit boundaries (sub-voxel). NaN if none.

u = u(:).'; u = u ./ max(norm(u), eps);
[M,N,P] = size(BW_metal);
xv = 1:N; yv = 1:M; zv = 1:P;

% sample t
t_samp = t_range(1):step_mm:t_range(2);
if numel(t_samp) < 2, t_samp = [t_range(1) t_range(2)]; end
Pk = p_start + t_samp(:).*u;

% RAS -> voxel (continuous, 1-based)
ijk_h = (Avox2ras \ [Pk, ones(size(Pk,1),1)].').';
v0 = ijk_h(:,1:3);
J = v0(:,1)+1; I = v0(:,2)+1; K = v0(:,3)+1;

% Trilinear sampling of metal mask
vals = interp3(xv, yv, zv, double(BW_metal), J, I, K, 'linear', 0);
inside = vals > 0.5;

% Find all inside runs; pick the **widest**
d = diff(inside);
rises = find(d == 1);
falls = find(d == -1);
if inside(1),  rises = [1; rises]; end
if inside(end), falls = [falls; numel(inside)-1]; end

if isempty(rises) || isempty(falls)
    enter_ras = [NaN NaN NaN]; exit_ras = enter_ras; return;
end

[~, pick] = max(t_samp(falls+1) - t_samp(rises));
ei = rises(pick); xo = falls(pick);

% Refine both boundaries with bisection on continuous mask value
t1 = refine_boundary(@(tt) interp_inside(BW_metal, Avox2ras, p_start + tt*u, M,N,P), t_samp(ei), t_samp(ei+1), 8);
t2 = refine_boundary(@(tt) interp_inside(BW_metal, Avox2ras, p_start + tt*u, M,N,P), t_samp(xo), t_samp(xo+1), 8);

enter_ras = p_start + t1*u;
exit_ras  = p_start + t2*u;
end

function inside = interp_inside(BW, Avox2ras, P_ras, M,N,P)
ijk_h = (Avox2ras \ [P_ras(:).', 1].').'; v0 = ijk_h(1:3);
J = v0(1)+1; I = v0(2)+1; K = v0(3)+1;
if J<1 || J>N || I<1 || I>M || K<1 || K>P, inside = false; return; end
inside = interp3(1:N,1:M,1:P,double(BW), J, I, K, 'linear', 0) > 0.5;
end

function t_star = refine_boundary(isInsideFn, ta, tb, nIters)
fa = isInsideFn(ta); fb = isInsideFn(tb);
if fa == fb, fb = isInsideFn(tb + 1e-3*max(1,abs(tb-ta))); end
for it=1:nIters
    tm = 0.5*(ta+tb);
    fm = isInsideFn(tm);
    if fa == fm, ta = tm; fa = fm; else, tb = tm; fb = fm; end
end
t_star = 0.5*(ta+tb);
end