function d2 = dist2_point_bezier(P, E, C, D, M)
[~,B] = bezier_samples(E,C,D,M);
d2 = min(pdist2(P, B, 'euclidean').^2, [], 2);
end