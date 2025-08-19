function d = pt_line_dist(P0, U, Q)
% Distance from points Q (N×3) to an infinite line P0 + t*U
% U must be unit length.

% enforce row vector and unit norm
U = U(:).';
U = U ./ max(norm(U), eps);

V = Q - P0;              % N×3
proj = (V * U.');        % N×1  (scalar projection)
par  = proj .* U;        % N×3  (component parallel to U)  <- implicit expansion
d    = vecnorm(V - par, 2, 2);
end