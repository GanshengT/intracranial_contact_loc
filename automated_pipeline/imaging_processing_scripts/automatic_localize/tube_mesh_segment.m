function [V,F] = tube_mesh_segment(E, D, r, nTheta, nSeg)
% Build a capped cylinder mesh around segment E->D
U = D - E; L = norm(U);
if L < eps, V = []; F = []; return; end
u = U / L;

% Orthonormal frame (v1,v2 orthogonal to u)
ref = [1 0 0]; if abs(dot(ref,u)) > 0.9, ref = [0 1 0]; end
v1 = cross(u, ref); v1 = v1 / max(norm(v1),eps);
v2 = cross(u, v1);  % already unit

t = linspace(0,1,nSeg+1);             % along axis
ang = linspace(0,2*pi,nTheta+1); ang(end) = [];  % around radius

% vertices along the side
V = zeros((nSeg+1)*nTheta + 2, 3);
for i = 1:numel(t)
    C = E + t(i)*U;
    for j = 1:numel(ang)
        idx = (i-1)*nTheta + j;
        V(idx,:) = C + r*(cos(ang(j))*v1 + sin(ang(j))*v2);
    end
end
% cap centers
idxCap1 = (nSeg+1)*nTheta + 1;  V(idxCap1,:) = E;
idxCap2 = (nSeg+1)*nTheta + 2;  V(idxCap2,:) = D;

% faces for side
F = [];
for i = 1:nSeg
    baseA = (i-1)*nTheta; baseB = i*nTheta;
    for j = 1:nTheta
        a1 = baseA + j;
        a2 = baseA + mod(j, nTheta) + 1;
        b1 = baseB + j;
        b2 = baseB + mod(j, nTheta) + 1;
        F = [F; a1 b1 b2; a1 b2 a2]; %#ok<AGROW> % two triangles per quad
    end
end
% faces for caps (fan)
for j = 1:nTheta
    j2 = mod(j, nTheta) + 1;
    % proximal cap
    F = [F; idxCap1, j2, j]; %#ok<AGROW>
    % distal cap ring indices
    dj  = nSeg*nTheta + j;
    dj2 = nSeg*nTheta + j2;
    F = [F; idxCap2, dj, dj2]; %#ok<AGROW>
end
end