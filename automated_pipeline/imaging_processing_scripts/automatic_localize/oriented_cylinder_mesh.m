function [V,F] = oriented_cylinder_mesh(C,u,len_mm,rad_mm,nTheta,nAx)
    u = u./max(norm(u),eps);
    % orthonormal basis (v1,v2) spanning plane ⟂ u
    v1 = null(u).';      % returns 2x3; take first row as v1, second as v2
    if size(v1,1) < 2
        V=[]; F=[]; return;
    end
    v2 = v1(2,:); v1 = v1(1,:);
    v1 = v1./norm(v1); v2 = v2./norm(v2);

    t = linspace(-len_mm/2, +len_mm/2, nAx+1);
    th = linspace(0, 2*pi, nTheta+1); th(end) = []; % avoid duplicate seam
    [T,TH] = ndgrid(t, th);
    R = rad_mm * ones(size(T));

    % points in RAS
    V = C + T(:).*u + R(:).*cos(TH(:)).*v1 + R(:).*sin(TH(:)).*v2;

    % faces (quads -> triangles)
    nZ = numel(t); nA = numel(th);
    toIdx = @(iz,ia) (iz-1)*nA + ia;
    F = [];
    for iz = 1:(nZ-1)
        for ia = 1:nA
            ia2 = mod(ia, nA) + 1;
            v11 = toIdx(iz,ia);
            v12 = toIdx(iz,ia2);
            v21 = toIdx(iz+1,ia);
            v22 = toIdx(iz+1,ia2);
            F = [F; v11 v21 v22; v11 v22 v12]; %#ok<AGROW>
        end
    end
end