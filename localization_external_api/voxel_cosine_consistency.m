function cosine_stats = voxel_cosine_consistency(P, percentile)
% For each voxel in P (Nx3), compute average pairwise cosine similarity 
% among the farthest "percentile" % of other voxels.
%
% Inputs:
%   P          - Nx3 point cloud
%   percentile - scalar in (0, 100), e.g., 50 means use farthest 50% of voxels
%
% Outputs:
%   cosine_stats.avg_cosine_similarity (Nx1)
%   cosine_stats.all_similarities (NxM) where M = num pairs used per voxel

assert(percentile > 0 && percentile < 100, 'percentile must be between 0 and 100');

N = size(P, 1);
avg_sim = zeros(N, 1);
all_sim = cell(N, 1);

for i = 1:N
    this_pt = P(i, :);
    others_idx = setdiff(1:N, i);
    others = P(others_idx, :);

    % Compute distances to all others
    D = sqrt(sum((others - this_pt).^2, 2));

    % Select top percentile (farthest points)
    num_select = max(3, round(percentile/100 * length(D)));
    [~, idx_sorted] = sort(D, 'descend');
    selected_idx = idx_sorted(1:num_select);

    selected_pts = others(selected_idx, :);
    vectors = bsxfun(@minus, selected_pts, this_pt); % Mx3

    % Normalize vectors
    norms = sqrt(sum(vectors.^2, 2)) + eps;
    Vnorm = bsxfun(@rdivide, vectors, norms); % Mx3 unit vectors

    % Compute pairwise cosine similarity among selected vectors
    cos_mat = Vnorm * Vnorm'; % MxM
    upper_vals = abs(cos_mat(triu(true(num_select), 1)));
    avg_sim(i) = mean(upper_vals);
    all_sim{i} = upper_vals;
end

cosine_stats.avg_cosine_similarity = avg_sim;
cosine_stats.all_similarities = all_sim;
end