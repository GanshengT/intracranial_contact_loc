function C = brain_centroid_ras(BW, Avox2ras0)
lin = find(BW);
if numel(lin) > 5e5, lin = lin(randperm(numel(lin), 5e5)); end
[i,j,k] = ind2sub(size(BW), lin);
V0 = [j i k] - 1;                          % [x y z] = [j i k] 0-based
P  = ([V0, ones(size(V0,1),1)] * Avox2ras0.');  % -> RAS
C  = mean(P(:,1:3),1);
end