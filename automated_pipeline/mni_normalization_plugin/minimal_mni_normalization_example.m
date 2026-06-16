%% Minimal example: normalize localized contacts to MNI or another reference brain
% This example uses a small BJH025 input dataset distributed with the plugin.
% Set environment variable RUN_MNI_REGISTRATION=true, or set run_registration
% below to true, to run ANTs registration and write transformed contacts.
% Leave it false to verify paths and inputs without launching ANTs.

clear;
clc;

plugin_dir = fileparts(mfilename('fullpath'));
addpath(genpath(plugin_dir));

subj_id = "BJH025";
subject_mri = fullfile(plugin_dir, 'example_data', subj_id, ...
    'imaging_process', 'SAG_T1_MPRAGE_converted_from_orig_mgz.nii');
contact_table_mat = fullfile(plugin_dir, 'example_data', subj_id, ...
    'imaging_process', 'contact_data_table_matched_bci2000_signal.mat');
output_dir = fullfile(plugin_dir, 'example_output', subj_id);

% This packaged NIfTI reference is copied from FreeSurfer:
% /Applications/freesurfer/8.1.0/average/mni_icbm152_nlin_asym_09c/
% Users can replace it with another reference anatomical image, such as a
% different MNI template, atlas-space brain, or study-specific reference MRI.
reference_mri = fullfile(plugin_dir, 'example_data', 'reference', ...
    'mni_icbm152_t1_tal_nlin_asym_09c.nii.gz');

% This packaged ANTs bin folder contains the runtime commands needed by the
% example. Replace it with your local ANTs installation if needed.
ants_bin_path = fullfile(plugin_dir, 'example_data', 'ants', 'bin');

run_registration = strcmpi(getenv('RUN_MNI_REGISTRATION'), 'true');

[contact_tbl, paths] = gt_normalize_contacts_to_reference_brain( ...
    subject_mri, contact_table_mat, reference_mri, output_dir, ...
    "SubjectID", subj_id, ...
    "AntsBinPath", ants_bin_path, ...
    "RunRegistration", run_registration);

disp(paths);
if run_registration
    disp(contact_tbl(:, {'ShankID', 'Reference_X', 'Reference_Y', 'Reference_Z'}));
else
    fprintf(['Dry run complete. Set run_registration = true to run ANTs ' ...
        'and write reference-space contact outputs.\n']);
end
