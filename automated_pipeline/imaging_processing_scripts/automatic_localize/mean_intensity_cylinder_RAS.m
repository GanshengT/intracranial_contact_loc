function mu = mean_intensity_cylinder_RAS(C,u,rad_mm,len_mm,Avox2ras0,ctVol,varargin)
% Robust contact score = trimmed mean(inside) - trimmed mean(annulus)
p = inputParser;
addParameter(p,'Trim',10);          % % trimmed at each tail
addParameter(p,'Annulus',[1.2 2.0]);% background ring radii (×rad)
addParameter(p,'Grid',[9 13 16]);   % [nr nz ntheta] samples
parse(p,varargin{:});
trim=p.Results.Trim; ann=p.Results.Annulus; g=p.Results.Grid;
nr=g(1); nz=g(2); ntheta=g(3);

% Orthonormal frame
u = u(:).'/max(norm(u),eps);
ref=[1 0 0]; if abs(dot(ref,u))>0.9, ref=[0 1 0]; end
v1 = cross(u,ref); v1 = v1/norm(v1+eps);
v2 = cross(u,v1);

% Sample points in RAS
r  = linspace(0,rad_mm,nr);
z  = linspace(-len_mm/2,+len_mm/2,nz);
th = linspace(0,2*pi,ntheta+1); th(end)=[];
[RR,ZZ,TT] = ndgrid(r,z,th);
Pin = C + ZZ(:).*u + RR(:).*(cos(TT(:)).*v1 + sin(TT(:)).*v2);

% background annulus
ra = linspace(ann(1)*rad_mm,ann(2)*rad_mm,max(3,ceil(nr/2)));
[RA,ZA,TA] = ndgrid(ra,z,th);
Pbg = C + ZA(:).*u + RA(:).*(cos(TA(:)).*v1 + sin(TA(:)).*v2);

% RAS -> voxel(1-based continuous) using 0-based Avox2ras
ijk_in = ras_to_ijk1(Avox2ras0, Pin);
ijk_bg = ras_to_ijk1(Avox2ras0, Pbg);

% Trilinear sampling
Iin = trilinear(double(ctVol), ijk_in);
Ibg = trilinear(double(ctVol), ijk_bg);

% Trimmed means
Iin = sort(Iin); Ibg = sort(Ibg);
tlo = floor(numel(Iin)*trim/200)+1; thi = numel(Iin)-floor(numel(Iin)*trim/200);
blo = floor(numel(Ibg)*trim/200)+1; bhi = numel(Ibg)-floor(numel(Ibg)*trim/200);
mu  = mean(Iin(tlo:thi),'omitnan') - mean(Ibg(blo:bhi),'omitnan');
end


% Trilinear interpolation in 1-based index space
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

% Tangent estimator from discrete centerline (linear interp)
function t = gradient_along_centerline(Cpts, s, s_query)
    % local finite diff on s
    s_query = max(s(1), min(s(end), s_query));
    x = interp1(s, Cpts(:,1), s_query, 'linear');
    y = interp1(s, Cpts(:,2), s_query, 'linear');
    z = interp1(s, Cpts(:,3), s_query, 'linear');
    % numerical gradient by small step (1 mm)
    h = 0.5;
    s0 = max(s(1), s_query - h);
    s1 = min(s(end), s_query + h);
    p0 = [interp1(s,Cpts(:,1),s0), interp1(s,Cpts(:,2),s0), interp1(s,Cpts(:,3),s0)];
    p1 = [interp1(s,Cpts(:,1),s1), interp1(s,Cpts(:,2),s1), interp1(s,Cpts(:,3),s1)];
    t = p1 - p0;
end

function [v1,v2,u] = axis_frame(u)
u = u(:).'; u = u / max(norm(u),eps);
ref = [1 0 0]; if abs(dot(ref,u))>0.9, ref=[0 1 0]; end
v1 = cross(u, ref); v1 = v1 / max(norm(v1),eps);
v2 = cross(u, v1);
end
