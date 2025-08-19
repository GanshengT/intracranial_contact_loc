function P_ras = ijk1_to_ras(Avox2ras, ijk1)
% ijk1: [N x 3] 1-based voxel indices
V0    = ijk1 - 1;
P_ras = ([V0, ones(size(V0,1),1)] * Avox2ras.'); % [N x 4]
P_ras = P_ras(:,1:3);
end