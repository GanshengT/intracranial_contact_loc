function visualize_pose_ct(C, u, rad_mm, len_mm, Avox2ras, ctVol, voxdims, thrHU, viz_half_mm, vec_len_mm, show_core_samples)
% Visualize candidate pose C (RAS), basis {u,v1,v2}, and nearby CT metal (>thrHU).
% C          : 1x3, RAS mm
% u          : 1x3, unit axis in RAS
% rad_mm     : contact radius (mm)
% len_mm     : contact length (mm)   (only used for sample illustration)
% Avox2ras   : 4x4 affine (voxel(0-based) -> RAS)
% ctVol      : MxNxP CT (HU)
% voxdims    : [dx dy dz] mm
% thrHU      : HU threshold to show metal blob (e.g., 3000)
% viz_half_mm: crop half-width around C in mm (e.g., 10–15)
% vec_len_mm : arrow length in mm (e.g., 8–12)
% show_core_samples : logical flag to scatter some core samples
    figure;

    if nargin < 11 || isempty(show_core_samples), show_core_samples = false; end

    C = C(:).'; 
    u = u(:).'; u = u./max(norm(u),eps);

    % Orthonormal basis in RAS
    [v1,v2] = orthobasis_row(u);

    % ----- crop a subvolume around C -----
    Ainv = inv(Avox2ras);
    ijk0   = [C,1]*Ainv.';     % 0-based voxel coords
    i0 = ijk0(1)+1; j0 = ijk0(2)+1; k0 = ijk0(3)+1;  % 1-based for indexing

    % convert mm half-width to voxels (each axis separately)
    half_i = max(2, round(viz_half_mm / max(voxdims(1),eps)));
    half_j = max(2, round(viz_half_mm / max(voxdims(2),eps)));
    half_k = max(2, round(viz_half_mm / max(voxdims(3),eps)));

    [m,n,p] = size(ctVol);
    Ir = max(1, floor(i0-half_i)) : min(m, ceil(i0+half_i));
    Jr = max(1, floor(j0-half_j)) : min(n, ceil(j0+half_j));
    Kr = max(1, floor(k0-half_k)) : min(p, ceil(k0+half_k));

    sub = ctVol(Ir, Jr, Kr);
    metal = sub > thrHU;

    % Build voxel-index grids for isosurface (meshgrid = [rows, cols, slices])
    [I,J,K] = meshgrid(Ir, Jr, Kr);    % note: meshgrid order is (Y,X,Z) visually; we pass I,J,K consistently
    % % isosurface on logical; take 0.5
    % try
    %     hs = patch(isosurface(I, J, K, permute(metal,[2 1 3]), 0.5));  % permute to keep (I=row, J=col) consistent
    % catch
    %     % If permute confuses, use simpler: build metal on the absolute grid directly
    %     hs = patch(isosurface(I, J, K, metal, 0.5));
    % end
    hs = patch(isosurface(I, J, K, metal, 0.5));
    set(hs,'EdgeColor','none','FaceAlpha',0.25,'FaceColor',[0.7 0.7 0.9]); % translucent lilac

    % Convert the vertices (voxel indices) -> RAS
    Vijk = hs.Vertices;                     % [i j k] (1-based)
    Vras = ([Vijk - 1, ones(size(Vijk,1),1)] * Avox2ras.').';
    Vras = Vras(1:3,:).';
    hs.Vertices = Vras;

    hold on; axis equal vis3d; grid on; box on;
    xlabel('X (RAS)'); ylabel('Y (RAS)'); zlabel('Z (RAS)');
    title('Pose debug: C, u, v1, v2 and local CT metal (> threshold)');

    % ----- draw C (center) -----
    plot3(C(1),C(2),C(3),'ko','MarkerFaceColor','y','MarkerSize',8,'DisplayName','C (center)');

    % ----- draw axis and basis vectors -----
    L = vec_len_mm;
    quiver3(C(1),C(2),C(3), L*u(1), L*u(2), L*u(3), 0, 'r','LineWidth',2, 'DisplayName','u (axis)');
    quiver3(C(1),C(2),C(3), L*v1(1),L*v1(2),L*v1(3),0, 'g','LineWidth',2, 'DisplayName','v1 (⊥)');
    quiver3(C(1),C(2),C(3), L*v2(1),L*v2(2),L*v2(3),0, 'b','LineWidth',2, 'DisplayName','v2 (⊥)');

    % ----- (optional) scatter a few core samples on the ring -----
    if show_core_samples
        th = linspace(0,2*pi,24);
        % a single axial slice near the center (t=0)
        r_core = rad_mm * 1.0;                  % show the physical core
        Pcore = C + r_core*cos(th(:))*v1 + r_core*sin(th(:))*v2;  % (24×3)
        scatter3(Pcore(:,1),Pcore(:,2),Pcore(:,3), 12, 'm', 'filled', 'DisplayName','core ring samples');
    end

    camlight headlight; lighting gouraud;
    legend('Location','best');
end