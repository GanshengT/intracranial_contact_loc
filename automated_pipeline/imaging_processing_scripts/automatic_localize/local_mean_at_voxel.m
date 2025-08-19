function m = local_mean_at_voxel(ctVol, i, j, k, rad_vox)
I = max(1, i-rad_vox(2)) : min(size(ctVol,1), i+rad_vox(2));
J = max(1, j-rad_vox(1)) : min(size(ctVol,2), j+rad_vox(1));
K = max(1, k-rad_vox(3)) : min(size(ctVol,3), k+rad_vox(3));
sub = double(ctVol(I,J,K));
m = mean(sub(:), 'omitnan');
end