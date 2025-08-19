function mask = cylinder_roi_mask(sz, Avox2ras, p0, p1, rad_mm)
% Exact cylinder (finite segment) around line p0->p1 in RAS (mm)
% sz: size(ctVol), Avox2ras: 4x4 vox->RAS, p0/p1: 1x3 or 3x1 (RAS mm), rad_mm: scalar radius
    p0 = p0(:).';  % force 1x3
    p1 = p1(:).';

    % Precompute voxel-center coordinates in RAS
    [I,J,K] = meshgrid(1:sz(1), 1:sz(2), 1:sz(3));     % voxel indices (1-based)
    V0  = [I(:)-1, J(:)-1, K(:)-1, ones(numel(I),1)];% 0-based for vox2ras1
    RAS = (V0 * Avox2ras.');                         % Nx4
    RAS = RAS(:,1:3);                                % Nx3 in mm

    % Vector math for distance to segment
    v   = p1 - p0;           % 1x3
    L2  = dot(v,v);          % squared length
    if L2 < eps
        % Degenerate: sphere around p0
        d = vecnorm(RAS - p0, 2, 2);
        mask = reshape(d <= rad_mm, sz);
        return;
    end
    w   = bsxfun(@minus, RAS, p0);   % Nx3
    t   = (w*v.')./L2;               % Nx1 projection param
    t   = max(0, min(1, t));         % clamp to segment
    proj = p0 + t.*v;                % Nx3 closest points on segment
    d    = vecnorm(RAS - proj, 2, 2);% Nx1 distances in mm

    mask = reshape(d <= rad_mm, sz);
end