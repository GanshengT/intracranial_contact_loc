function [ii,jj,kk, dmin] = k_nearest_to_point(mask, Xras, P0, K)
    idx = find(mask);
    if isempty(idx), ii=[]; jj=[]; kk=[]; dmin=[]; return; end
    [im,jm,km] = ind2sub(size(mask), idx);
    % coords of candidates
    Pc = [ squeeze(Xras(sub2ind([size(mask) 3], im,jm,km, ones(size(im))*1))), ...
           squeeze(Xras(sub2ind([size(mask) 3], im,jm,km, ones(size(im))*2))), ...
           squeeze(Xras(sub2ind([size(mask) 3], im,jm,km, ones(size(im))*3))) ];
    % distances to P0 (seed RAS)
    d = sqrt(sum((Pc - P0(:)'.^1).^2, 2));
    [~,ord] = sort(d,'ascend');
    ord = ord(1:min(K,numel(ord)));
    ii = im(ord); jj = jm(ord); kk = km(ord);
    dmin = d(ord);
end