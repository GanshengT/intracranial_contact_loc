function d2 = dist2_point_line_segment(P, E, D)
% squared distance from each row of P to line segment E--D
U = D - E; L2 = sum(U.^2);
t = ((P - E) * U.') / max(L2, eps);
t = min(max(t,0),1);
Proj = E + t.*U;
d2 = sum((P - Proj).^2, 2);
end