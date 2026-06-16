%% External API example:
% This may be useful for researcher working with previous data normalized
% at MNI space
% call the local MNI/reference normalization plugin
% The canonical implementation lives in
%  automated_pipeline/mni_normalization_plugin.

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
plugin_dir = fullfile(repo_root, 'automated_pipeline', 'mni_normalization_plugin');
addpath(genpath(plugin_dir));

subj_id = "BJH025";
subject_mri = fullfile(plugin_dir, 'example_data', subj_id, ...
    'imaging_process', 'SAG_T1_MPRAGE_converted_from_orig_mgz.nii');
contact_table_mat = fullfile(plugin_dir, 'example_data', subj_id, ...
    'imaging_process', 'contact_data_table_matched_bci2000_signal.mat');
output_dir = fullfile(plugin_dir, 'example_output', subj_id);

% This packaged example reference is copied from FreeSurfer. Replace it
% with any other reference anatomical image appropriate for your target.
reference_mri = fullfile(plugin_dir, 'example_data', 'reference', ...
    'mni_icbm152_t1_tal_nlin_asym_09c.nii.gz');

% This packaged ANTs bin folder contains the runtime commands needed by the
% example. Replace it with your local ANTs installation if needed.
ants_bin_path = fullfile(plugin_dir, 'example_data', 'ants', 'bin');
% run_registration = strcmpi(getenv('RUN_MNI_REGISTRATION'), 'true');
run_registration = true;
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
