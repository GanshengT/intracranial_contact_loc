function P = sample_line(p0, p1, step_mm)
% uniform samples between p0 and p1 (RAS, mm)
    L   = norm(p1 - p0);
    n   = max(2, ceil(L/step_mm));
    t   = linspace(0,1,n).';
    P   = p0 + t.*(p1 - p0);
end