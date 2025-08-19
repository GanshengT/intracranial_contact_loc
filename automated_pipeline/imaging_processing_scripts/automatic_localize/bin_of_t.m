function tb = bin_of_t(t, edges)
    tb = max(1, min(numel(edges)-1, 1 + floor((t - edges(1)) / (edges(2)-edges(1)))));
end
