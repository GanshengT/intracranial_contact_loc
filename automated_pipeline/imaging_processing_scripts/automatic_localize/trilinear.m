function val = trilinear(V, ijk1)
% V(i,j,k), ijk1(:,1)=i,(:,2)=j,(:,3)=k  (all 1-based, continuous)
sz = size(V); i = ijk1(:,1); j = ijk1(:,2); k = ijk1(:,3);
i = max(1,min(sz(1), i)); j = max(1,min(sz(2), j)); k = max(1,min(sz(3), k));
i0=floor(i); i1=min(i0+1,sz(1));  di=i-i0;
j0=floor(j); j1=min(j0+1,sz(2));  dj=j-j0;
k0=floor(k); k1=min(k0+1,sz(3));  dk=k-k0;
idx = @(ii,jj,kk) sub2ind(sz,ii,jj,kk);
c000=V(idx(i0,j0,k0)); c100=V(idx(i1,j0,k0)); c010=V(idx(i0,j1,k0)); c110=V(idx(i1,j1,k0));
c001=V(idx(i0,j0,k1)); c101=V(idx(i1,j0,k1)); c011=V(idx(i0,j1,k1)); c111=V(idx(i1,j1,k1));
c00 = c000.*(1-di) + c100.*di;  c01 = c001.*(1-di) + c101.*di;
c10 = c010.*(1-di) + c110.*di;  c11 = c011.*(1-di) + c111.*di;
c0  = c00.*(1-dj) + c10.*dj;    c1  = c01.*(1-dj) + c11.*dj;
val = c0.*(1-dk) + c1.*dk;
end