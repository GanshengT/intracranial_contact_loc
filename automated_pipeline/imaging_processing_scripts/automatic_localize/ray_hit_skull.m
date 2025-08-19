function entry_ras = ray_hit_skull(BW_skull, Avox2ras, P0, U, step_mm)
    % March from P0 along +U until the first skull voxel; returns [NaN NaN NaN] if none.
    U = U(:).'; U = U ./ max(norm(U), eps);
    sz = size(BW_skull);
    
    % conservative bound (200 mm)
    nSteps = ceil(200/step_mm);
    entry_ras = [NaN NaN NaN];
    for k = 0:nSteps
        Pk  = P0 + (k*step_mm).*U;
        ijk = ras2ijk1(Avox2ras, Pk);
        if any(ijk<1) || ijk(1)>sz(1) || ijk(2)>sz(2) || ijk(3)>sz(3), continue; end
        if BW_skull(ijk(1),ijk(2),ijk(3))
            entry_ras = Pk; return;
        end
    end
end
