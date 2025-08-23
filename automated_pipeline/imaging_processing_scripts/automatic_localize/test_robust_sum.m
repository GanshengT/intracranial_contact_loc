%TEST_ROBUST_SUM Unit tests for robust_sum function.
%   Verifies handling of quantile bounds and invalid values.

assert(isequal(robust_sum([1 2 3 4], 0.5), 7));
assert(isequal(robust_sum([1 2 3 4], 2), 10));
assert(isequal(robust_sum([1 2 3 4], -1), 4));
assert(isequal(robust_sum([1 2 NaN Inf 5], 0.5), 7));
