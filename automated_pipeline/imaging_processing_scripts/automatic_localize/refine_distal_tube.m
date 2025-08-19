function [distal_refined, best_score] = refine_distal_tube(P, entry_ras, e1, distal_init, r_mm, maxGrow)
    % 1-D search extending distal by up to maxGrow mm to maximize tube overlap.
    L0 = norm(distal_init - entry_ras);
    if L0 < eps, distal_refined = distal_init; best_score = tube_overlap_score(P, entry_ras, distal_init, r_mm); return; end
    
    cand = linspace(L0, L0 + maxGrow, 11); % 0..maxGrow in 1-mm-ish steps
    best_score = -inf; distal_refined = distal_init;
    for L = cand
        dpt = entry_ras + L * e1;
        sc  = tube_overlap_score(P, entry_ras, dpt, r_mm);
        if sc > best_score
            best_score = sc; distal_refined = dpt;
        end
    end
end