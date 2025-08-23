function v = robust_sum(vals, q)
    %ROBUST_SUM Sum the largest fraction of values while ignoring invalid entries.
    %   V = ROBUST_SUM(VALS, Q) returns the sum of the largest Q fraction of
    %   the elements in VALS, where Q is between 0 and 1. Values that are Inf
    %   or NaN are ignored. Q values outside the [0,1] range are clamped so
    %   the function never requests more elements than exist, avoiding errors
    %   from MAXK.

    vals = vals(:);
    vals = vals(~isinf(vals) & ~isnan(vals));
    if isempty(vals), v = -Inf; return; end
    n = numel(vals);
    k = min(n, max(1, ceil(q * n)));
    v = sum(maxk(vals, k));
end
