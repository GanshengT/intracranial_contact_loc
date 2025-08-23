function r = safe_ratio(a,b)
    %SAFE_RATIO Compute ratio while avoiding division by zero.
    %   R = SAFE_RATIO(A,B) returns A divided by B, but if B is zero the
    %   result is set to zero instead of inf or NaN.  Negative denominators
    %   are handled normally.

    if b == 0
        r = 0;
    else
        r = a / b;
    end
end
