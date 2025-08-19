function show_three_axial_slices_before_after(C_ref_ras, ct_before, Avox2ras_before, ...
                                               ct_after,  Avox2ras_after, figtitle)
    arguments
        C_ref_ras (1,3) double
        ct_before  double
        Avox2ras_before (4,4) double
        ct_after   double
        Avox2ras_after  (4,4) double
        figtitle   char = 'Axial slices near shank (before vs after)'
    end

    % Helper: RAS -> voxel index (0-based), then to MATLAB 1-based
    ras2vox = @(A, P) ([P 1] / A.').';   % returns homogeneous column vector
    uvw0_b  = ras2vox(Avox2ras_before, C_ref_ras);  uvw0_b = uvw0_b(1:3);
    uvw0_a  = ras2vox(Avox2ras_after , C_ref_ras);  uvw0_a = uvw0_a(1:3);

    k0_b = round(uvw0_b(3)) + 1;   % axial index in BEFORE
    k0_a = round(uvw0_a(3)) + 1;   % axial index in AFTER

    % Clamp to valid range and build triplets
    kz_b = max(2, min(size(ct_before,3)-1, k0_b)) + (-1:1);
    kz_a = max(2, min(size(ct_after ,3)-1, k0_a)) + (-1:1);

    clim = [0 3000];

    figure('Color','w','Name',figtitle); 
    tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for i=1:3
        nexttile(tl, i);
        imagesc(squeeze(ct_before(:,:,kz_b(i)))'); axis image off; colormap(gray); caxis(clim);
        title(sprintf('Before  k=%d', kz_b(i)));
    end
    for i=1:3
        nexttile(tl, 3+i);
        imagesc(squeeze(ct_after(:,:,kz_a(i)))'); axis image off; colormap(gray); caxis(clim);
        title(sprintf('After   k=%d', kz_a(i)));
    end
    title(tl, figtitle);
end