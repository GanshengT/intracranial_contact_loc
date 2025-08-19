function [B1, B2, t1, t2] = brain_boundary_hits(BW, Avox2ras0, p0, u, t_range, step_mm)

u = u(:).'; u = u ./ max(norm(u), eps);
[M,N,P] = size(BW);
xv = 1:N; yv = 1:M; zv = 1:P;

% Sample along line
t_samp = t_range(1):step_mm:t_range(2);
if numel(t_samp) < 2, t_samp = [t_range(1) t_range(2)]; end
Pk = p0 + t_samp(:).*u;

% RAS -> voxel (continuous, 1-based)
ijk_h = (Avox2ras0 \ [Pk, ones(size(Pk,1),1)].').';
v0 = ijk_h(:,1:3);
J = v0(:,1)+1; I = v0(:,2)+1; K = v0(:,3)+1;

% Trilinear sampling for robust inside test
vals = interp3(xv, yv, zv, double(BW), J, I, K, 'linear', 0);
inBrain = vals > 0.5;

% Find widest inside run 
d = diff(inBrain);
rises = find(d == 1);
falls = find(d == -1);
if inBrain(1), rises = [1; rises]; end
if inBrain(end), falls = [falls; numel(inBrain)-1]; end
if isempty(rises) || isempty(falls), B1=nan(1,3); B2=B1; t1=nan; t2=nan; return; end

% choose pair (rise, fall) with max span
[~,pick] = max(t_samp(falls+1) - t_samp(rises));
enter_idx = rises(pick);
exit_idx  = falls(pick);

% refine with bisection
t1 = refine_boundary(@(tt) interp_inside(BW,Avox2ras0,p0 + tt*u,M,N,P), t_samp(enter_idx), t_samp(enter_idx+1), 8);
t2 = refine_boundary(@(tt) interp_inside(BW,Avox2ras0,p0 + tt*u,M,N,P), t_samp(exit_idx),  t_samp(exit_idx+1),  8);
B1 = p0 + t1*u;
B2 = p0 + t2*u;
end

function inside = interp_inside(BW, Avox2ras0, P_ras, M,N,P)
ijk_h = (Avox2ras0 \ [P_ras(:).', 1].').'; v0 = ijk_h(1:3);
J = v0(1) + 1; I = v0(2) + 1; K = v0(3) + 1;
if J<1 || J>N || I<1 || I>M || K<1 || K>P
    inside = false; return;
end
inside = interp3(1:N,1:M,1:P,double(BW), J, I, K, 'linear', 0) > 0.5;
end

function t_star = refine_boundary(isInsideFn, ta, tb, nIters)
fa = isInsideFn(ta); fb = isInsideFn(tb);
if fa == fb
    eps_t = 1e-3 * max(1, abs(tb - ta));
    fb = isInsideFn(tb + eps_t);
end
for it=1:nIters
    tm = 0.5*(ta+tb);
    fm = isInsideFn(tm);
    if fa == fm
        ta = tm; fa = fm;
    else
        tb = tm; fb = fm;
    end
end
t_star = 0.5*(ta+tb);
end

