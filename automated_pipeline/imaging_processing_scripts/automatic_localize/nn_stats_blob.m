function [muD, sigD, stats, d_vec] = nn_stats_blob(G, Xras, varargin)
% NN_STATS_BLOB
% Compute mean/SD of inter-voxel distances for the current blob.
% Modes:
%   'pairwise' : mean/sd over ALL pairwise distances (O(N^2)); trims optional.
%   'kNN'      : mean/sd over the K nearest distances per voxel (robust & scalable).
%
% Outputs:
%   muD, sigD : mean and std of the distance distribution used
%   stats     : struct with fields: Mode, N, K, TrimPct, UsedN, PairwiseComputed (bool)

    ip = inputParser; ip.CaseSensitive=false;
    addParameter(ip,'Mode','kNN',@(s)ischar(s)&&ismember(lower(s),{'pairwise','knn'}));
    addParameter(ip,'K',8,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'TrimPct',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<50);   % per tail
    addParameter(ip,'MaxPairwiseN',3000,@(x)isnumeric(x)&&isscalar(x)&&x>=100); % cap for O(N^2)
    addParameter(ip,'SampleN',2000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);     % when subsampling
    parse(ip,varargin{:});
    Mode         = lower(ip.Results.Mode);
    K            = ip.Results.K;
    TrimPct      = ip.Results.TrimPct;
    MaxPairwiseN = ip.Results.MaxPairwiseN;
    SampleN      = ip.Results.SampleN;

    idx = find(G);
    n = numel(idx);
    if n < 2, muD=0; sigD=0; stats=struct('Mode',Mode,'N',n,'K',K,'TrimPct',TrimPct,'UsedN',n,'PairwiseComputed',false); return; end

    [ii,jj,kk] = ind2sub(size(G), idx);
    P = [ squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*1))), ...
          squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*2))), ...
          squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*3))) ];

    switch Mode
        case 'pairwise'
            if n <= MaxPairwiseN
                d = pdist(P, 'euclidean');      % condensed vector, length n*(n-1)/2
                usedN = n; pairwiseOK = true;
            else
                % Subsample to control O(N^2) cost
                rp = randperm(n, min(SampleN, n));
                d = pdist(P(rp,:), 'euclidean');
                usedN = numel(rp); pairwiseOK = false;
            end

            % Optional trimming
            if TrimPct > 0
                d = sort(d);
                L = numel(d);
                iL = floor(TrimPct/100 * L) + 1;
                iU = ceil((1 - TrimPct/100) * L);
                d = d(max(1,iL):min(L,iU));
            end
            d_vec = d; 
            muD  = mean(d);
            sigD = std(d);

            stats = struct('Mode','pairwise','N',n,'K',NaN,'TrimPct',TrimPct, ...
                           'UsedN',usedN,'PairwiseComputed',pairwiseOK);

        case 'knn'
            % Compute kNN distances per point (exclude self)
            % For large n, this is still heavy, but much cheaper than full pairwise.
            D = pdist2(P, P);
            D(1:n+1:end) = inf;                % remove self
            Kuse = min(K, n-1);
            Dsort = mink(D, Kuse, 2);          % n x Kuse (requires R2018b+)

            % Optional trimming across the K neighbors
            if TrimPct > 0 && Kuse > 2
                Dsort = sort(Dsort, 2);
                tL = floor(TrimPct/100 * Kuse) + 1;
                tU = ceil((1 - TrimPct/100) * Kuse);
                Dtrim = Dsort(:, max(1,tL):min(Kuse,tU));
            else
                Dtrim = Dsort;
            end

            vec = Dtrim(:);
            d_vec = Dtrim(:); 
            muD  = mean(vec);
            sigD = std(vec);

            stats = struct('Mode','kNN','N',n,'K',Kuse,'TrimPct',TrimPct, ...
                           'UsedN',n,'PairwiseComputed',false);
    end
end