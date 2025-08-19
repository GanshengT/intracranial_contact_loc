function ijk1 = ras_to_ijk1(Avox2ras0, rasNx3)
% Avox2ras maps [x0 y0 z0 1]^T (0-based) -> RAS. Invert, then 0->1 based.
% Returns [i j k] (row, col, slice) as 1-based floating indices for MATLAB.
R = Avox2ras0;
xyz0 = (R \ [rasNx3, ones(size(rasNx3,1),1)].').';  % -> [x0 y0 z0]
j0 = xyz0(:,1);  i0 = xyz0(:,2);  k0 = xyz0(:,3);   % x->col(j), y->row(i), z->slice(k)
ijk1 = [i0+1, j0+1, k0+1];
end