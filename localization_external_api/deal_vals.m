function y = deal_vals(ss, F, a, b)
    mask = (ss>=a & ss<=b);
    y    = -Inf(size(ss));
    if any(mask), y(mask) = F(ss(mask)); end
end