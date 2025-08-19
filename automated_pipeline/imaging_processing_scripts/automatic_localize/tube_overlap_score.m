function s = tube_overlap_score(P, entry_ras, distal_ras, r_mm)
    % Fraction of blob points within a tube of radius r_mm around the line entry->distal
    U = distal_ras - entry_ras; U = U ./ max(norm(U), eps);
    V = P - entry_ras;                   % Nx3
    proj = (V * U.');                    % Nx1
    par  = proj .* U;                    % Nx3
    d    = vecnorm(V - par, 2, 2);       % Nx1 perp distances
    % Also clip longitudinally to segment [0, L]
    L = max(norm(distal_ras - entry_ras), eps);
    in_seg = proj >= 0 & proj <= L;
    s = mean((d <= r_mm) & in_seg);
end
