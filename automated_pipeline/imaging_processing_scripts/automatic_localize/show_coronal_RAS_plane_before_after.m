function show_coronal_RAS_plane_before_after(C_ref_ras, ct_before, Avox2ras_before, ...
                                             ct_after,  Avox2ras_after, figtitle)
    arguments
        C_ref_ras (1,3) double
        ct_before  double
        Avox2ras_before (4,4) double
        ct_after   double
        Avox2ras_after  (4,4) double
        figtitle   char = 'Coronal plane in RAS (before vs after)'
    end

    % RAS plane settings (in mm)
    width_mm  = 250;   % R/L span
    height_mm = 250;   % S/I span
    step_mm   = 0.5;   % sampling step

    % Build a coronal plane in RAS: span X (R/L) and Z (S/I); normal is Y (A/P)
    u = -width_mm/2  : step_mm :  width_mm/2;   % along R/L (X)
    v = -height_mm/2 : step_mm :  height_mm/2;  % along S/I (Z)
    [U,V] = meshgrid(u,v);

    Xras = C_ref_ras(1) + U;              % vary in X
    Yras = C_ref_ras(2) + 0*U;            % fixed coronal plane (Y const)
    Zras = C_ref_ras(3) + V;              % vary in Z

    % Helper: RAS -> voxel (column output)
    ras2vox = @(A, P) ([P 1] / A.').';  % P is 1x3, returns 4x1; take 1:3

    % Vectorize: map the plane to voxel coords for each volume
    function [I,J,K] = rasgrid2vox(A, X, Y, Z)
        n = numel(X);
        XYZ1 = [X(:) Y(:) Z(:) ones(n,1)];
        uvw  = (XYZ1 / A.').';      % 4 x n
        uvw  = uvw(1:3,:).';        % n x 3
        I = reshape(uvw(:,1)+1, size(X));   % 1-based for MATLAB
        J = reshape(uvw(:,2)+1, size(Y));
        K = reshape(uvw(:,3)+1, size(Z));
    end

    [Ib,Jb,Kb] = rasgrid2vox(Avox2ras_before, Xras, Yras, Zras);
    [Ia,Ja,Ka] = rasgrid2vox(Avox2ras_after , Xras, Yras, Zras);

    % Sample both volumes on the same RAS plane
    % Use linear interp, extrapolate to NaN
    before_slice = interp3(ct_before, Jb, Ib, Kb, 'linear', NaN);
    after_slice  = interp3(ct_after , Ja, Ia, Ka, 'linear', NaN);

    clim = [0 3000];

    figure('Color','w','Name',figtitle);
    t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(u, v, before_slice); axis image tight; colormap(gray); caxis(clim);
    xlabel('R/L (mm)'); ylabel('S/I (mm)'); title('Before (coronal RAS)');

    nexttile;
    imagesc(u, v, after_slice);  axis image tight; colormap(gray); caxis(clim);
    xlabel('R/L (mm)'); ylabel('S/I (mm)'); title('After (coronal RAS)');

    title(t, figtitle);
end