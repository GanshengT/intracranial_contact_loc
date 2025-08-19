function [edges, cap] = summarize_fill_metal(cand_base, t_field, TBinMM)
    t = t_field(cand_base);
    if isempty(t), edges=[0 1]; cap=0; return; end
    t0 = min(t); t1 = max(t);
    nb = max(1, ceil((t1 - t0)/TBinMM));
    edges = linspace(t0, t1, nb+1);
    cap = histcounts(t, edges);  % how much METAL is available per slab
end