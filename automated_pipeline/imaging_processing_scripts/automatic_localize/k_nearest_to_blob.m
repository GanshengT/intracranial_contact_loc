function [ic,jc,kc, dmin_vec, davg_vec, dmed_vec] = k_nearest_to_blob( ...
        G, eligible, Xras, KNearest, varargin)
% K_NEAREST_TO_BLOB
% Returns K candidates with the smallest MIN distance to the blob,
% plus per-candidate average and median distances to (subsampled) blob points.
%
% Name-Value:
%   'SampleN'   (2000) : subsample this many blob points to bound memory/time
%   'AvgMode'    ('mean') : 'mean' | 'median' | 'trimmed'
%   'TrimPct'   (10)   : if AvgMode='trimmed', percentage to trim on each tail

    ip = inputParser; ip.CaseSensitive=false;
    addParameter(ip,'SampleN',2000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
    addParameter(ip,'AvgMode','mean',@(s)ischar(s) && ismember(lower(s),{'mean','median','trimmed'}));
    addParameter(ip,'TrimPct',10,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<50);
    parse(ip,varargin{:});
    SampleN = ip.Results.SampleN;
    AvgMode = lower(ip.Results.AvgMode);
    TrimPct = ip.Results.TrimPct;

    idxE = find(eligible);
    if isempty(idxE)
        ic=[]; jc=[]; kc=[]; dmin_vec=[]; davg_vec=[]; dmed_vec=[];
        return;
    end

    % Candidate RAS points
    [ie,je,ke] = ind2sub(size(G), idxE);
    Pc = [ squeeze(Xras(sub2ind([size(G) 3], ie,je,ke, ones(size(ie))*1))), ...
           squeeze(Xras(sub2ind([size(G) 3], ie,je,ke, ones(size(ie))*2))), ...
           squeeze(Xras(sub2ind([size(G) 3], ie,je,ke, ones(size(ie))*3))) ];

    % Blob RAS points (possibly subsampled for speed)
    [ig,jg,kg] = ind2sub(size(G), find(G));
    if isempty(ig)
        dmin_vec = inf(size(idxE));
        davg_vec = inf(size(idxE));
        dmed_vec = inf(size(idxE));
    else
        Pg_all = [ squeeze(Xras(sub2ind([size(G) 3], ig,jg,kg, ones(size(ig))*1))), ...
                   squeeze(Xras(sub2ind([size(G) 3], ig,jg,kg, ones(size(ig))*2))), ...
                   squeeze(Xras(sub2ind([size(G) 3], ig,jg,kg, ones(size(ig))*3))) ];
        if size(Pg_all,1) > SampleN
            rp = randperm(size(Pg_all,1), SampleN);
            Pg = Pg_all(rp, :);
        else
            Pg = Pg_all;
        end

        % Pairwise distances (bounded by SampleN)
        D = pdist2(Pc, Pg);   % size: (#cand) x (#blob_sample)

        % Per-candidate summaries
        dmin_vec = min(D, [], 2);

        switch AvgMode
            case 'mean'
                davg_vec = mean(D, 2);
            case 'median'
                davg_vec = median(D, 2);
            case 'trimmed'
                qL = TrimPct/100; qU = 1 - qL;
                % sort each row once; then trim
                Dsort = sort(D, 2);
                n = size(Dsort,2);
                iL = max(1, floor(qL*n));
                iU = min(n, ceil(qU*n));
                davg_vec = mean(Dsort(:, iL:iU), 2);
            otherwise
                davg_vec = mean(D, 2);
        end

        dmed_vec = median(D, 2);  % handy if you want it
    end

    % Rank by min-distance (stable with avg as tiebreaker if you want)
    [~, ord] = sort(dmin_vec, 'ascend');
    K = min(KNearest, numel(ord));
    ord = ord(1:K);

    idxE     = idxE(ord);
    dmin_vec = dmin_vec(ord);
    davg_vec = davg_vec(ord);
    dmed_vec = dmed_vec(ord);

    [ic,jc,kc] = ind2sub(size(G), idxE);
end