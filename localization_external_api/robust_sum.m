function v = robust_sum(vals, q)
    vals = vals(:);
    vals = vals(~isinf(vals) & ~isnan(vals));
    if isempty(vals), v = -Inf; return; end
    k = max(1, ceil(q*numel(vals)));
    v = sum(maxk(vals, k));
end