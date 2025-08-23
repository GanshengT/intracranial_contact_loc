%TEST_SAFE_RATIO Unit tests for safe_ratio function.
%   Exercises basic cases including zero and negative denominators.

assert(isequal(safe_ratio(6, 3), 2));
assert(isequal(safe_ratio(6, -3), -2));
assert(isequal(safe_ratio(6, 0), 0));
