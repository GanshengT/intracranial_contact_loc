function ctr = blob_centroid_ras(G, Xras)
    idx = find(G);
    if isempty(idx), ctr = squeeze(Xras(1,1,1,:)).'; return; end
    [ii,jj,kk] = ind2sub(size(G), idx);
    P = [ squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*1))), ...
          squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*2))), ...
          squeeze(Xras(sub2ind([size(G) 3], ii,jj,kk, ones(size(ii))*3))) ];
    ctr = mean(P,1);
end