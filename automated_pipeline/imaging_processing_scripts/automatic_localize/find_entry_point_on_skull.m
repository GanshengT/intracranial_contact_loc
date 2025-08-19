function entry_ras = first_skull_before_t(BW_skull, Avox2ras, P0, U, t_stop, step_mm)
    sz = size(BW_skull);
    nSteps = max(1, ceil(t_stop/step_mm));
    for k = 0:nSteps
        Pk = P0 + (k*step_mm).*U;
        ijk = ras2ijk1(Avox2ras, Pk);
        if any(ijk<1) || ijk(1)>sz(1) || ijk(2)>sz(2) || ijk(3)>sz(3), continue; end
        if BW_skull(ijk(1),ijk(2),ijk(3))
            entry_ras = Pk; return;
        end
    end
    % fallback: try along full line
    nSteps = ceil(200/step_mm); % 200 mm safety
    for k = 0:nSteps
        Pk = P0 + (k*step_mm).*U;
        ijk = ras2ijk1(Avox2ras, Pk);
        if any(ijk<1) || ijk(1)>sz(1) || ijk(2)>sz(2) || ijk(3)>sz(3), continue; end
        if BW_skull(ijk(1),ijk(2),ijk(3))
            entry_ras = Pk; return;
        end
    end
    warning('No skull intersection found; using P0 as entry.'); entry_ras = P0;
end