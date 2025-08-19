function tube_mask = tube_mask_from_params(ctr, dstar, t0, t1, r_mm, Xras, Domain)
% TUBE_MASK_FROM_PARAMS  Voxelize a tube defined by (ctr, dstar, [t0,t1], r_mm).
% If Domain (logical) is provided, compute only on that subset for speed.

    if nargin < 7 || isempty(Domain)
        Domain = true(size(Xras,1), size(Xras,2), size(Xras,3));
    end
    idx = find(Domain);
    [ii,jj,kk] = ind2sub(size(Domain), idx);

    d = dstar(:).'/max(norm(dstar),eps);

    % coords relative to ctr for Domain voxels
    X = [ squeeze(Xras(sub2ind([size(Domain) 3], ii,jj,kk, ones(size(ii))*1))) - ctr(1), ...
          squeeze(Xras(sub2ind([size(Domain) 3], ii,jj,kk, ones(size(ii))*2))) - ctr(2), ...
          squeeze(Xras(sub2ind([size(Domain) 3], ii,jj,kk, ones(size(ii))*3))) - ctr(3) ];

    tloc = X*d(:);
    cx   = cross(X, repmat(d,numel(idx),1), 2);
    rloc = sqrt(sum(cx.^2,2));

    ok = (tloc >= t0) & (tloc <= t1) & (rloc <= r_mm);
    tube_mask = false(size(Domain));
    tube_mask(idx(ok)) = true;
end