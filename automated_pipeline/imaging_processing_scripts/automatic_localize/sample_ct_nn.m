function HU = sample_ct_nn(vol, Avox2ras, Xras)
% Nearest-voxel HU at RAS points Xras (N×3), robust to axis-order issues.
    Ainv = inv(Avox2ras);
    uvw0 = [Xras, ones(size(Xras,1),1)] * Ainv.';   % -> [i j k] in 0-based
    Iq = uvw0(:,1) + 1;                              % rows
    Jq = uvw0(:,2) + 1;                              % cols
    Kq = uvw0(:,3) + 1;                              % slices
    % interp take ijk
    HU = interp3(double(vol), Iq, Jq, Kq, 'nearest', NaN);  % X=J, Y=I, Z=K !!
end