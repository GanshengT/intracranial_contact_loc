function tube_mask = voxelized_tube(ctr, d, r_mm, TPadMM, G, Xras, t_field)
    % Build a tube around axis (ctr,d), radius r_mm, length spanning G +/- TPad
    if any(G(:))
        tG = t_field(G); t0 = min(tG)-TPadMM; t1 = max(tG)+TPadMM;
    else
        t_all = t_field(:); t0 = min(t_all); t1 = max(t_all);
    end
    % project voxel centers onto (ctr,d)
    % shift origin to ctr
    X = Xras; X(:,:,:,1) = X(:,:,:,1) - ctr(1);
    X(:,:,:,2) = X(:,:,:,2) - ctr(2);
    X(:,:,:,3) = X(:,:,:,3) - ctr(3);
    % longitudinal and radial
    tloc = X(:,:,:,1)*d(1) + X(:,:,:,2)*d(2) + X(:,:,:,3)*d(3);
    crossX = cat(4, ...
        X(:,:,:,2)*d(3) - X(:,:,:,3)*d(2), ...
        X(:,:,:,3)*d(1) - X(:,:,:,1)*d(3), ...
        X(:,:,:,1)*d(2) - X(:,:,:,2)*d(1));
    rloc = sqrt(sum(crossX.^2,4));
    tube_mask = (tloc >= t0) & (tloc <= t1) & (rloc <= r_mm);
end
