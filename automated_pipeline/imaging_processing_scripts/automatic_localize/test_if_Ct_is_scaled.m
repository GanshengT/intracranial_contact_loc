% load two cts, get transformation
% internal control, not related to the main function% 
ct1_filename = '/Users/ganshengtan/Library/CloudStorage/Box-Box/BJH079/IMAGING/NIfTI/CT_BONE/CT_BONE.img';
orginal_ct_filename = '/Users/ganshengtan/Library/CloudStorage/Box-Box/BJH079/IMAGING/DICOM/CT_BONE/IM-0001-0001-0001.dcm';

spm('defaults','fmri'); spm_jobman('initcfg');

% to nifti
[dcm_folder,~,~] = fileparts(orginal_ct_filename);
dcm_list = dir(fullfile(dcm_folder, '*.dcm'));
if isempty(dcm_list), error('No DICOMs found in %s', dcm_folder); end

hdr = spm_dicom_headers(fullfile({dcm_list.folder}, {dcm_list.name}));
out = spm_dicom_convert(hdr, 'all', 'flat', 'nii');   % returns struct with .files
src_files = out.files;                                % cellstr of produced .nii
src_filename = src_files{1};                          % take first volume of the series

fprintf('DICOM converted -> NIfTI: %s\n', src_filename);

% coreg
Vref = spm_vol(ct1_filename);
Vsrc_before = spm_vol(src_filename);
Mref        = Vref.mat;           % voxel -> world (mm)
Msrc_before = Vsrc_before.mat;

fixed  = spm_read_vols(Vref);
moving = spm_read_vols(Vsrc_before);

% Normalize intensities a bit (helps MI)
fixed  = mat2gray(fixed(isfinite(fixed)));   fixed  = reshape(fixed,  Vref.dim);
moving = mat2gray(moving(isfinite(moving))); moving = reshape(moving, Vsrc_before.dim);

% affine transform
% Mattes MI (default metric), tune optimizer for 3D
[optimizer, metric] = imregconfig('multimodal');
optimizer.GrowthFactor     = 1.05;
optimizer.Epsilon          = 1.5e-6;
optimizer.InitialRadius    = 6.25e-3;
optimizer.MaximumIterations= 300;

tform_vox = imregtform(moving, fixed, 'affine', optimizer, metric, ...
    'PyramidLevels', 4);  % maps moving (src vox) -> fixed (ref vox)

T_vox = tform_vox.T;   % 4x4 in voxel index coordinates

%update header ----------
Delta_affine_mm = Mref * T_vox * inv(Msrc_before);   % world(mm): old_src_world -> ref_world

Msrc_after = Delta_affine_mm * Msrc_before;          % new header mat for src
spm_get_space(Vsrc_before.fname, Msrc_after);        

% Report 12-parameter decomposition ----------
p = spm_imatrix(Delta_affine_mm);  % [Tx Ty Tz Rx Ry Rz Zx Zy Zz Kx Ky Kz]
tx=p(1); ty=p(2); tz=p(3);
rx=p(4); ry=p(5); rz=p(6);
zx=p(7); zy=p(8); zz=p(9);
kx=p(10); ky=p(11); kz=p(12);

fprintf('\n=== Affine coreg (src -> ref) ===\n');
disp('Δ_affine (world mm) = Mref * T_vox * inv(Msrc_before):');
disp(Delta_affine_mm);
fprintf('Translation (mm): [%.3f %.3f %.3f]\n', tx,ty,tz);
fprintf('Rotation (rad) : [%.6f %.6f %.6f]\n', rx,ry,rz);
fprintf('Rotation (deg) : [%.3f %.3f %.3f]\n', [rx ry rz]*180/pi);
fprintf('Scales        : [%.6f %.6f %.6f]\n', zx,zy,zz);
fprintf('Shears        : [%.6f %.6f %.6f]\n', kx,ky,kz);

% decompose: [Tx Ty Tz Rx Ry Rz Zx Zy Zz Kx Ky Kz]
% get transform parameter
% delta diagonal not 0 reflects results from rotation 

