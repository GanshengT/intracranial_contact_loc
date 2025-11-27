% main function - minimal working example
% this demonstrate what minimum input files is needed
% full functinality requires using output from GT's coregistration package

im_struct = load('BJH075_imaging_processing_info_struct.mat');
subj_id = 'BJH075';
im.ct_volume_registered_to_mri = im_struct.ct_volume_registered_to_mri;
% this is the nifti file, you can do load nifti for the coregistered Ct
im.coregistered_to_mri_trajectory = im_struct.coregistered_to_mri_trajectory;
% this is traj defined by two points, you can get from ROSA files, can
% manually defined (not recommend because neurologists have different
% naming scheme). In GT's coregistration, this is adjusted from relative
% locations between ROSA's executed trajs, not disclosed here
im.electrode_info = im_struct.electrode_info;
% read from OR sheet, reading function for standard WashU OR sheet is available:
% g.tan@wustl.edu

%   im        : struct with fields (all in coregistered-CT→MRI space)
% ct_volume_registered_to_mri can be obtained using MRIread(CT.nii),
% MRIread from freesurfer/Matlab
%                 .ct_volume_registered_to_mri.vol       [I×J×K] CT (numeric)
%                 .ct_volume_registered_to_mri.vox2ras1  4×4 (voxel1→RAS)
%                 .ct_volume_registered_to_mri.vox2ras0  4×4 (voxel0→RAS)
%                 .coregistered_to_mri_trajectory        table/struct with:
%                       .trajectory_id{n} (char/string)
%                       .trajectory(n).start, .trajectory(n).end (1×3 RAS mm)
%                 .electrode_info table with columns:
%                       n_contact_str, n_contact, etc. (row-aligned to trajectories)
% this is important, because we require on manufacture spec, if electrodes
% are from different manufactures. Please use a list such as ['PMT',
% 'adtech', ...], currently support adtech, PMT and DIXI, you can add more
electrode_manufacture = 'DIXI';
% 5 either from Gt's sugmentation or from freesurfer (not freesurfer might
% call it brainmask.mgz), I used automatic segmentation
bm_path = 'BJH075_brainmask.auto.mgz';
% you can set the folder path for storing the results here, for example,
% curr_dir = pwd;
% then input 'curr_dir', curr_dir to the function.

[contact_tbl, shank_model_all, paths] = localize_electrode_gt_small_example(im, subj_id, "electrode_manufacture", electrode_manufacture, ...
    "bm_path", bm_path);
