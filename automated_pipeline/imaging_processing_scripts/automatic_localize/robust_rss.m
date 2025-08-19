function RSS = robust_rss(d2, trimPct)
d2 = sort(d2(:));
n  = numel(d2);
k0 = floor(n*trimPct/100)+1;
k1 = n - floor(n*trimPct/100);
RSS = sum(d2(max(1,k0):max(k0,k1)));
end