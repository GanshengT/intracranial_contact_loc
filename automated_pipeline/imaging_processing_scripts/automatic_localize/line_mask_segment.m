function [B1, B2, t1, t2] = line_mask_segment(BW, Avox2ras, p0, u, t_range, step_mm)
% Find the two boundary intersections of line p0 + t u with a binary mask BW (CT grid).
% Returns RAS boundary points B1 (enter) and B2 (exit) with sub-voxel refinement.

u = u(:).'; u = u ./ max(norm(u), eps);
t_samp = t_range(1):step_mm:t_range(2);
inBrain = false(size(t_samp));
sz = size(BW);

for k = 1:numel(t_samp)
    Pk = p0 + t_samp(k)*u;
    ijk = ras2ijk1(Avox2ras, Pk);
    if any(ijk<1) || ijk(1)>sz(1) || ijk(2)>sz(2) || ijk(3)>sz(3)
        inBrain(k) = false;
    else
        inBrain(k) = BW(ijk(1),ijk(2),ijk(3));
    end
end

% transitions false->true (enter) and true->false (exit)
enter_idx = find(~inBrain(1:end-1) & inBrain(2:end), 1, 'first');
exit_idx  = find(inBrain(1:end-1) & ~inBrain(2:end), 1, 'last');

if isempty(enter_idx) || isempty(exit_idx)
    B1 = [NaN NaN NaN]; B2 = [NaN NaN NaN]; t1 = NaN; t2 = NaN; return;
end

% refine each boundary with binary search on t
t1 = refine_boundary(@(tt) mask_at(BW, Avox2ras, p0 + tt*u, sz), t_samp(enter_idx), t_samp(enter_idx+1), 6);
t2 = refine_boundary(@(tt) mask_at(BW, Avox2ras, p0 + tt*u, sz), t_samp(exit_idx),  t_samp(exit_idx+1),  6);

B1 = p0 + t1*u; 
B2 = p0 + t2*u;
end