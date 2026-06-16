function [contact_tbl, paths] = gt_normalize_contacts_to_reference_brain( ...
    subject_mri, contact_table_input, reference_mri, output_dir, opts)
%GT_NORMALIZE_CONTACTS_TO_REFERENCE_BRAIN Normalize contact coordinates to a reference brain.
%
% This plugin registers a subject anatomical MRI to a reference brain image
% with ANTs and applies the validated affine-only coordinate transform to
% localized intracranial contacts. The reference image can be FreeSurfer's
% MNI152 template, another MNI image, an atlas-space brain, or a study
% specific reference MRI.
%
% Example:
%   [contact_tbl, paths] = gt_normalize_contacts_to_reference_brain( ...
%       subject_mri, contact_table_mat, reference_mri, output_dir, ...
%       "SubjectID", "BJH025", "AntsBinPath", ants_bin_path);
%
% Required contact columns:
%   ShankID, X_aligned_to_brightest_voxel, Y_aligned_to_brightest_voxel,
%   Z_aligned_to_brightest_voxel

arguments
    subject_mri (1,1) string
    contact_table_input
    reference_mri (1,1) string
    output_dir (1,1) string
    opts.SubjectID (1,1) string = ""
    opts.AntsBinPath (1,1) string = ""
    opts.OverwriteExisting (1,1) logical = false
    opts.RunRegistration (1,1) logical = true
    opts.TransformMode (1,1) string {mustBeMember(opts.TransformMode, "affine_only")} = "affine_only"
    opts.XColumn (1,1) string = "X_aligned_to_brightest_voxel"
    opts.YColumn (1,1) string = "Y_aligned_to_brightest_voxel"
    opts.ZColumn (1,1) string = "Z_aligned_to_brightest_voxel"
    opts.WriteDatFiles (1,1) logical = true
    opts.Verbose (1,1) logical = true
end

subject_mri = char(subject_mri);
reference_mri = char(reference_mri);
output_dir = char(output_dir);

assert(exist(subject_mri, 'file') == 2, ...
    'Subject MRI does not exist: %s', subject_mri);
assert(exist(reference_mri, 'file') == 2, ...
    'Reference MRI does not exist: %s', reference_mri);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

if strlength(opts.SubjectID) == 0
    [~, inferred_id] = fileparts(fileparts(fileparts(subject_mri)));
    opts.SubjectID = string(inferred_id);
end

paths = struct();
paths.output_dir = output_dir;
paths.subject_mri = subject_mri;
paths.reference_mri = reference_mri;
paths.output_prefix = fullfile(output_dir, 'reference_');
paths.affine_file = fullfile(output_dir, 'reference_0GenericAffine.mat');
paths.warped_image = fullfile(output_dir, 'reference_Warped.nii.gz');
paths.inverse_warped_image = fullfile(output_dir, 'reference_InverseWarped.nii.gz');
paths.contact_patient_csv = fullfile(output_dir, 'contacts_patient_space.csv');
paths.contact_reference_csv = fullfile(output_dir, 'contacts_reference_affine_only.csv');
paths.dat_folder = fullfile(output_dir, 'reference_contact_coordinate');

contact_tbl = loadContactTable(contact_table_input);
contact_tbl = validateContactTable(contact_tbl, opts);

if opts.RunRegistration
    ants_bin_path = char(opts.AntsBinPath);
    if isempty(ants_bin_path)
        ants_bin_path = findAntsBinPath();
    end
    setupAntsPath(ants_bin_path);

    if opts.OverwriteExisting || exist(paths.affine_file, 'file') ~= 2 || ...
            exist(paths.warped_image, 'file') ~= 2
        % linear registration for efficiency (non linear > 2hrs)
        cmd = sprintf('%s -d 3 -t a -f %s -m %s -o %s', ...
            shellQuote(fullfile(ants_bin_path, 'antsRegistrationSyN.sh')), ...
            shellQuote(reference_mri), shellQuote(subject_mri), ...
            shellQuote(paths.output_prefix));
        if opts.Verbose
            fprintf('Running ANTs registration:\n%s\n', cmd);
        end
        [status, result] = system(cmd);
        if opts.Verbose
            fprintf('%s\n', result);
        end
        assert(status == 0, 'ANTs registration failed.');
    elseif opts.Verbose
        fprintf('Reusing existing registration in %s.\n', output_dir);
    end
else
    if opts.Verbose
        fprintf(['RunRegistration is false. Path and contact-table checks ' ...
            'completed without ANTs registration.\n']);
    end
    return
end

writetable(contact_tbl, paths.contact_patient_csv);

assert(exist(paths.affine_file, 'file') == 2, ...
    'Missing affine transform: %s', paths.affine_file);
assert(exist(paths.warped_image, 'file') == 2, ...
    'Missing warped image: %s', paths.warped_image);

affine_struct = load(paths.affine_file);
affine_matrix = ea_antsmat2mat( ...
    affine_struct.AffineTransform_double_3_3, affine_struct.fixed);

coords = [ ...
    contact_tbl.(opts.XColumn), ...
    contact_tbl.(opts.YColumn), ...
    contact_tbl.(opts.ZColumn), ...
    ones(height(contact_tbl), 1)];
transformed_coords = coords * affine_matrix';

contact_tbl.Reference_X = transformed_coords(:, 1);
contact_tbl.Reference_Y = transformed_coords(:, 2);
contact_tbl.Reference_Z = transformed_coords(:, 3);
contact_tbl.TransformMode = repmat(opts.TransformMode, height(contact_tbl), 1);
contact_tbl.SubjectID = repmat(opts.SubjectID, height(contact_tbl), 1);

% Future nonlinear point-transform option:
% The validated default is affine_only. To apply the full ANTs affine + warp
% to contact points, replace the multiplication above with
% antsApplyTransformsToPoints. ANTs point transforms are applied in the
% opposite direction of image transforms. For subject-space points to a
% reference image, write physical coordinates, handle RAS/LPS consistently,
% and use [reference_0GenericAffine.mat,1] with reference_1InverseWarp.nii.gz.

writetable(contact_tbl, paths.contact_reference_csv);
if opts.WriteDatFiles
    writeContactDatFiles(contact_tbl, paths.dat_folder);
end
validateCoordinateRange(contact_tbl, opts.SubjectID);

if opts.Verbose
    fprintf('Saved reference-space contacts: %s\n', paths.contact_reference_csv);
end
end

function contact_tbl = loadContactTable(contact_table_input)
if istable(contact_table_input)
    contact_tbl = contact_table_input;
    return
end

contact_table_path = char(string(contact_table_input));
assert(exist(contact_table_path, 'file') == 2, ...
    'Contact table input does not exist: %s', contact_table_path);
loaded = load(contact_table_path);

if isfield(loaded, 'contact_info_seeg_matched_bci2000_signal')
    contact_tbl = loaded.contact_info_seeg_matched_bci2000_signal;
elseif isfield(loaded, 'contact_data_table')
    contact_tbl = loaded.contact_data_table;
else
    fields = strjoin(fieldnames(loaded), ', ');
    error(['Could not find contact_info_seeg_matched_bci2000_signal or ' ...
        'contact_data_table in %s. Found fields: %s'], contact_table_path, fields);
end
end

function contact_tbl = validateContactTable(contact_tbl, opts)
required_vars = ["ShankID", opts.XColumn, opts.YColumn, opts.ZColumn];
missing_vars = setdiff(required_vars, string(contact_tbl.Properties.VariableNames));
assert(isempty(missing_vars), ...
    'Contact table is missing required columns: %s', strjoin(missing_vars, ', '));

for coord_var = [opts.XColumn, opts.YColumn, opts.ZColumn]
    values = contact_tbl.(coord_var);
    assert(isnumeric(values), 'Coordinate column %s must be numeric.', coord_var);
    assert(all(isfinite(values)), ...
        'Coordinate column %s contains NaN or Inf.', coord_var);
end
end

function setupAntsPath(ants_bin_path)
assert(exist(ants_bin_path, 'dir') == 7, ...
    'ANTs bin folder does not exist: %s', ants_bin_path);
setenv('PATH', [ants_bin_path ':' getenv('PATH')]);

required_bins = ["antsRegistrationSyN.sh", "antsRegistration", ...
    "PrintHeader", "antsApplyTransforms", "antsApplyTransformsToPoints"];
for i_bin = 1:numel(required_bins)
    executable_path = fullfile(ants_bin_path, char(required_bins(i_bin)));
    assert(exist(executable_path, 'file') == 2, ...
        'Missing ANTs executable: %s', executable_path);
    [status, result] = system(sprintf('chmod +x %s', shellQuote(executable_path)));
    assert(status == 0, 'Could not chmod +x %s: %s', executable_path, result);
end
end

function ants_bin_path = findAntsBinPath()
[status, result] = system('command -v antsRegistrationSyN.sh');
assert(status == 0 && ~isempty(strtrim(result)), ...
    ['Could not find antsRegistrationSyN.sh on PATH. ' ...
     'Pass "AntsBinPath", "/path/to/ants/bin" to the function.']);
ants_bin_path = fileparts(strtrim(result));
end

function writeContactDatFiles(contact_tbl, dat_folder)
if ~exist(dat_folder, 'dir')
    mkdir(dat_folder);
end

unique_shank_ids = unique(string(contact_tbl.ShankID), 'stable');
for i_shank = 1:numel(unique_shank_ids)
    shank_id = unique_shank_ids(i_shank);
    filtered_tbl = contact_tbl(string(contact_tbl.ShankID) == shank_id, :);
    dat_path = fullfile(dat_folder, ...
        sprintf('%s_reference_affine.dat', sanitizeFileName(shank_id)));

    fid = fopen(dat_path, 'w');
    assert(fid ~= -1, 'Could not open %s for writing.', dat_path);
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '\n');
    for i_contact = 1:height(filtered_tbl)
        fprintf(fid, '%f %f %f\n', ...
            filtered_tbl.Reference_X(i_contact), ...
            filtered_tbl.Reference_Y(i_contact), ...
            filtered_tbl.Reference_Z(i_contact));
    end
    fprintf(fid, 'info\n');
    fprintf(fid, 'numpoints %d\n', height(filtered_tbl));
    fprintf(fid, 'useRealRAS 1\n');
    clear cleanup;
end
end

function validateCoordinateRange(contact_tbl, subj_id)
reference_xyz = [ ...
    contact_tbl.Reference_X, ...
    contact_tbl.Reference_Y, ...
    contact_tbl.Reference_Z];
fprintf('\nReference-space coordinate sanity check for %s\n', subj_id);
fprintf('Reference X range: %.2f to %.2f\n', ...
    min(reference_xyz(:,1)), max(reference_xyz(:,1)));
fprintf('Reference Y range: %.2f to %.2f\n', ...
    min(reference_xyz(:,2)), max(reference_xyz(:,2)));
fprintf('Reference Z range: %.2f to %.2f\n', ...
    min(reference_xyz(:,3)), max(reference_xyz(:,3)));
assert(all(isfinite(reference_xyz(:))), ...
    'Reference-space contact coordinates contain NaN or Inf for %s.', subj_id);
if any(abs(reference_xyz(:)) > 130)
    warning(['Reference coordinates for %s exceed +/-130 mm. ' ...
        'Check transform direction, coordinate space, and reference image.'], subj_id);
end
end

function mat = ea_antsmat2mat(afftransform, m_Center)
mat = [reshape(afftransform(1:9), [3, 3])', afftransform(10:12)];
m_Translation = mat(:, 4);
mat = [mat; [0, 0, 0, 1]];

for i = 1:3
    m_Offset(i) = m_Translation(i) + m_Center(i); %#ok<AGROW>
    for j = 1:3
        m_Offset(i) = m_Offset(i) - (mat(i, j) * m_Center(j));
    end
end

mat(1:3, 4) = m_Offset;
mat = inv(mat);
mat = mat .* [ ...
     1  1 -1 -1; ...
     1  1 -1 -1; ...
    -1 -1  1  1; ...
     1  1  1  1];
end

function safe_name = sanitizeFileName(name)
safe_name = char(name);
safe_name = regexprep(safe_name, '[/:*?"<>|]', '_');
safe_name = strtrim(safe_name);
if isempty(safe_name)
    safe_name = 'unnamed_shank';
end
end

function quoted = shellQuote(path_name)
quoted = ['''' strrep(char(path_name), '''', '''"''"''') ''''];
end
