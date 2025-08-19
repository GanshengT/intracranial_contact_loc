function t = t_param(P0, U, Q)
% P0: 1x3, U: 1x3 (unit), Q: Nx3
U = U(:).';                         % ensure row vector
U = U ./ max(norm(U), eps);         % unit
V = Q - P0;                         % Nx3
t = V * U.';                        % Nx1 scalar projection (no dot/broadcasting)
end