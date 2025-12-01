function [contact_tbl, shank_model_all] = gt_automatic_localize_electrode(imaging_processing_info_struct, subj_id, opts)
% This script is a packaged function for automatic localization
% This script outputs the contact coordinates and label
% Figures showing the contact location will be generated
% Update - Aug 17, 2025
% update note, we will employ even-interval constraint
% we will use linear prior, than we add one flex pt, see if differ than
% prior, if so, add one flex pt to increase nonlinearity.

% In the current pipeline, we are doing electrode localization on
% coregistered CT, the coregistration is rigid transformation, in this case
% electrode geometry is preserved. Affine transformation can be considered
% if we are doing longitudinal studies where brain deforms

% In general, we should avoid warping the electrode, to preserve geometry,
% we can warp MRi into electrode space. Or use dots represent rigid
% transformation of the electrode to get brain parcellation label. Also, 
% I recommend against reslicing ct for electrode localization.
%
% Inputs
%   im        : struct with fields (all in coregistered-CT→MRI space)
%                 .ct_volume_registered_to_mri.vol       [I×J×K] CT (numeric)
%                 .ct_volume_registered_to_mri.vox2ras1  4×4 (voxel1→RAS)
%                 .ct_volume_registered_to_mri.vox2ras0  4×4 (voxel0→RAS)
%                 .coregistered_to_mri_trajectory        table/struct with:
%                       .trajectory_id{n} (char/string)
%                       .trajectory(n).start, .trajectory(n).end (1×3 RAS mm)
%                 .electrode_info table with columns:
%                       n_contact_str, n_contact, etc. (row-aligned to trajectories)
%   subj_id   : char/string, used in filenames
%   opts      : struct of options (all optional)
%                 .manufacturer         default 'PMT', can be DIXI or
%                 adtech
%                 .report_dir           default [pwd filesep subj_id filesep 'imaging_process']
%                 .temp_dir             default same as report_dir
%                 .use_shank_model      default true
%                 .use_predefined_loc   default false (expects DATs in opts.manual_dir)
%                 .manual_dir           path with manually_identified_macro_contacts/*.dat
%                 .roi_radius_mm        default 4
%                 .metal_thr_override   scalar (use directly if given)
%                 .skull_thr_override   scalar (use directly if given)
%                 .save_intermediate    default true (pdf/mp4 in report_dir)
%                 .verbose              default true
%                 .manufacturer_map     containers.Map or struct: per-shank manufacturer override
%                 .subject_overrides    function handle @(subj_id, traj_id, start,end,par) -> struct of overrides
%                                        (e.g., to adjust ROI radius, thresholds, or endpoints)
%                 .seed_params          struct to override seed picking and growth params
%
% Outputs
%   contact_tbl        : table of contact centers and attributes (one row/contact)
%   shank_model_all    : struct array of shank models and per-contact results
%
% Notes
% Disabling any micro wire adjustment as this function is for macro
% localization
% Example
%   [T, S] = gt_automatic_localize_electrodes(imStruct, 'BJH033');
% bm_path is always required, if no bm_path is provided, a reference brain
% mask will be used, which is not optimal!

arguments
    imaging_processing_info_struct
    subj_id
    opts.curr_dir (1,1) string = pwd   % <--- NEW: default = pwd

    % System
    opts.freesurfer_matlab_path   (1,1) string  = "/Applications/freesurfer/8.1.0/matlab"
    opts.freesurfer_fsfast_path   (1,1) string  = "/Applications/freesurfer/8.1.0/fsfast/toolbox"
    opts.bm_path (1,1) string = 'reference_brainmask.auto.mgz';
    opts.use_individual_brainmask (1,1) logical = true

    % Electrode configuration
    opts.electrode_manufacture string = "PMT"
    opts.UseShankModel            (1,1) logical = true
    opts.UsePredefinedLoc         (1,1) logical = false
    opts.ManualDir                (1,1) string  = ""

    % ----- Geometry parameters (flat names, user-overridable) -----
    opts.adtech_BF_diameter        (1,1) double = 1.28
    opts.PMT_diameter              (1,1) double = 0.8
    opts.PMT_contact_length        (1,1) double = 2
    opts.adtech_BF_contact_length  (1,1) double = 1.57
    opts.DIXI_am_diameter          (1,1) double = 0.8
    opts.DIXI_am_contact_length    (1,1) double = 2
    opts.DIXI_am_spacing           (1,1) double = 3.5
    opts.DIXI_cm_larger_spacing    (1,1) double = 9
    opts.DIXI_bm_larger_spacing    (1,1) double = 13

    % ----- Electrode spacing models  -----
    opts.elec_model_dixi_15cm   (:,1) double = [0 3.5 7 10.5 14 27 30.5 34 37.5 41 54 57.5 61 64.5 68].'
    opts.elec_model_dixi_18cm   (:,1) double = [0 3.5 7 10.5 14 17.5 30.5 34 37.5 41 44.5 48 61 64.5 68 71.5 75 78.5].'
    opts.elec_model_dixi_am     (:,1) double = (0:3.5:59.5).'
    opts.elec_model_adtech      (:,1) double = (0:5:35).'
    opts.elec_model_dixi_mme    (:,1) double = [0 4 8 32 36 40].'
    opts.sub_vol_range_mm       (1,1) double = 1.57

    % Spatial and runtime options, we disable these functions in this
    % wrapper, so wont matter, keep them here for next iterations
    opts.RoiRadiusMM              (1,1) double  = 4
    opts.ReportDir                (1,1) string  = ""
    opts.SeedParams               struct       = struct()
    opts.Verbose                  (1,1) logical = true
end

addpath(opts.freesurfer_matlab_path)
% this is for faster
addpath(opts.freesurfer_fsfast_path)
% electrode_manufacture = opts.electrode_manufacture ;
n_elec = height(imaging_processing_info_struct.electrode_info);

if isscalar(opts.electrode_manufacture)
    % replicate the single manufacturer for all electrodes
    electrode_manufacture = repmat(opts.electrode_manufacture, n_elec, 1);
elseif isvector(opts.electrode_manufacture)
    % Case 2: user supplied list — must match electrode count
    if numel(opts.electrode_manufacture) ~= n_elec
        error("Mismatch: number of manufacturers (%d) ≠ number of electrodes (%d).", ...
              numel(opts.electrode_manufacture), n_elec);
    else
        electrode_manufacture = opts.electrode_manufacture(:);
    end
else
    error("opts.electrode_manufacture must be a string or string array.");
end


opts.electrode_manufacture = electrode_manufacture;

% if opts.Verbose
%     fprintf('Electrode manufacturer(s): %s (applied to %d electrodes)\n', ...
%         opts.electrode_manufacture(1), n_elec);
% end

% get parameters either default or user-defined
adtech_BF_diameter = opts.adtech_BF_diameter;
PMT_diameter = opts.PMT_diameter;
PMT_contact_length = opts.PMT_contact_length; % mm
adtech_BF_contact_length = opts.adtech_BF_contact_length; % mm
DIXI_am_diameter = opts.DIXI_am_diameter;
DIXI_am_contact_length = opts.DIXI_am_contact_length; % mm
DIXI_am_spacing = opts.DIXI_am_spacing;
DIXI_cm_larger_spacing = opts.DIXI_cm_larger_spacing;
DIXI_bm_larger_spacing = opts.DIXI_bm_larger_spacing;
use_predefined_location = opts.UsePredefinedLoc;
use_shank_model = opts.UseShankModel;
% research shank name - disabled
% DIXI_research_shank_id = {"B'", "G'"};
% we will consider the 3 tetrodes located at the first macro contact
% no research electrode
% research_shank_id = {};
% manufacture spect, might be useful
elec_model_dixi_15cm = opts.elec_model_dixi_15cm;
elec_model_dixi_18cm = opts.elec_model_dixi_18cm;
elec_model_dixi_am = opts.elec_model_dixi_am;
elec_model_adtech = opts.elec_model_adtech;
elec_model_dixi_mme = opts.elec_model_dixi_mme;

% disabled, different wrapper for localizing research shank - micro wires
% research_shank_id = {'A''^MTG-sts-Amy', 'B^MTG-aHc'};
% make the string array into char array
% research_shank_id = cellfun(@char, research_shank_id, 'UniformOutput', false);

sub_vol_range_mm = opts.sub_vol_range_mm; 
% mm unit, we want the subvol engulge the contact
% interp_ratio = 10; 
% % this parameter is no longer used
% to adjust macro contact location to the brightest voxel (most likely
% contact location), we will break a voxel into 10*10*10 small voxels.
% the current method to get objective macro electrode location:
% gradually decrease thresholding intensity until we found a connected
% volume
% continuous increase the volume if the connected volume does not meet the
% minimum_volume_for_contact
% also disabled this functionality for faster performance
% intensity_search_resolution = 1;
% macro_contact_diameter = 1.28;
% macro_contact_length = 1.57;
% macro_contact_volume_mm_cube = macro_contact_diameter^2 * pi * macro_contact_length;
% micro_contact_from_first_macro = 4; 
% minimum_voxel_for_micro_contact = 2; % this is empirical, adjust if needed
% we assume micro contact is 3mm from the first macro
% isovalue_to_get_surface = 0.5;
% % 0.5 is chosen for getting surface for binary volume
% % we will use nearest neighbor for label assignment, we first define the
% % volume for neighbor to be considered
% neighbor_distance_requirement = sqrt(macro_contact_diameter^2 + ...
%     macro_contact_length^2); % mm
% shank_start_from_last_macro_contact = 10; %mm, for plotting shank purpose
% shank_tip_diameter = 0.2;
% shank_tip_length = 1;
% the tip is a cylinder of 0.2mm diameter of 1mm length
% micro_contact_id = 99;

% default set at argument
% curr_dir = pwd;
curr_dir = opts.curr_dir;
temp_folder = fullfile(curr_dir, 'imaging_processing_temp_folder');

% Create it if it doesn't exist, for storing cache, again no need for this
% wrapper, only helpful for debugging
if ~exist(temp_folder, 'dir')
    mkdir(temp_folder);
end


automated_aligned_to_brightest_voxel_macro_contacts_folder = fullfile(curr_dir, ...
    '/automated_aligned_to_brightest_voxel_macro_contacts');
if ~exist(automated_aligned_to_brightest_voxel_macro_contacts_folder, 'dir')
    mkdir(automated_aligned_to_brightest_voxel_macro_contacts_folder);
end

automated_shank_center_line_folder = fullfile(curr_dir, ...
    '/automated_shank_center_lines');
if ~exist(automated_shank_center_line_folder, 'dir')
    mkdir(automated_shank_center_line_folder);
end

% save contacts' location visualization - disabled here
% report_contacts_related_surf_folder = fullfile(curr_dir, ...
%     '/report_contacts_related_surf_folder');
% if ~exist(report_contacts_related_surf_folder, 'dir')
%     mkdir(report_contacts_related_surf_folder);
% end

% region of interest - disabled, this is for GT's parcellation, one can use
% freesurfer as alternatives, freesurfer also provides some standard
% parcellation atlas
% we will focus on aparc - Desikan/Killiany parcellation
% region_of_interest = {'Left-Cerebral-White-Matter', 'Left-Thalamus', ...
%     'Left-Caudate', 'Left-Putamen', ...
%     'Left-Pallidum', 'Left-Hippocampus', 'Left-Amygdala', 'Left-Accumbens-area', ...
%     'Left-VentralDC', 'ctx-lh-bankssts', 'ctx-lh-caudalanteriorcingulate', ...
%     'ctx-lh-caudalmiddlefrontal', 'ctx-lh-corpuscallosum', 'ctx-lh-cuneus', ...
%     'ctx-lh-entorhinal', 'ctx-lh-fusiform', 'ctx-lh-inferiorparietal', ...
%     'ctx-lh-inferiortemporal', 'ctx-lh-isthmuscingulate', ...
%     'ctx-lh-lateraloccipital', 'ctx-lh-lateralorbitofrontal', ...
%     'ctx-lh-lingual', 'ctx-lh-medialorbitofrontal', 'ctx-lh-middletemporal', ...
%     'ctx-lh-parahippocampal', 'ctx-lh-paracentral', 'ctx-lh-parsopercularis', ...
%     'ctx-lh-parsorbitalis', 'ctx-lh-parstriangularis', 'ctx-lh-pericalcarine', ...
%     'ctx-lh-postcentral', 'ctx-lh-posteriorcingulate', 'ctx-lh-precentral', ...
%     'ctx-lh-precuneus', 'ctx-lh-rostralanteriorcingulate', ...
%     'ctx-lh-rostralmiddlefrontal', 'ctx-lh-superiorfrontal', ...
%     'ctx-lh-superiorparietal', 'ctx-lh-superiortemporal', ...
%     'ctx-lh-supramarginal', 'ctx-lh-frontalpole', 'ctx-lh-temporalpole', ...
%     'ctx-lh-transversetemporal', 'ctx-lh-insula'};
% region_of_interest_hippo_subsegment = {...
%     'presubiculum-head', 'subiculum-head', 'CA1-head', 'CA3-head', ...
%     'CA4-head', 'GC-ML-DG-head', 'molecular_layer_HP-head', 'Lateral-nucleus', ...
%     'Basal-nucleus', 'Central-nucleus', 'Medial-nucleus', 'Cortical-nucleus'...
%     'Accessory-Basal-nucleus'
%     };


% get CT intensity threshold, for bone and electrode volume
ctVol = imaging_processing_info_struct.ct_volume_registered_to_mri.vol;      % 3D numeric
Avox2ras1 = imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras1; % 4x4 affine (voxel->RAS)
Avox2ras0 = imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras0; % 4x4 affine (voxel->RAS)

v = double(ctVol(:));
v = v(isfinite(v));
v = v(v > 0);                      % drop air if desired
% reduced peak detection, 2100 is metal in most Cts
thr=2100;

% viz to confirm - deleted


% viz 3d construction of the bone + electrode vol
BW = ctVol >= thr;
% (Optional) remove tiny specks to de-noise the surface
% Use a size that’s small relative to electrode/bone blobs; tune as needed.
% because electrodes are separated, we should not remove small objects
% BW = bwareaopen(BW, 500);  % removes 3D components with <500 voxels, removve small objects ~ noises
% optional, could potential better represent the realistic shank
% BW_metal = imclose(BW_metal, strel('sphere', 1));     % tidy edges (optional)


% Build voxel index grid (1-based)
[m,n,p] = size(BW);
[X, Y, Z] = meshgrid(1:n, 1:m, 1:p);

% no need for viz,
% figure('Color','w'); hold on;
% h = patch(isosurface(X, Y, Z, BW, 0.5));             % use meshgrid coords
% isonormals(X, Y, Z, double(ctVol), h);               % normals from intensity
% set(h, 'EdgeColor','none', 'FaceColor',[0.75 0.78 0.85], 'FaceAlpha',0.25);



% get skull mask
ctNoMetal = double(ctVol);
ctNoMetal(ctVol >= thr) = NaN; 
w = ctNoMetal(:);
w = w(isfinite(w) & w > 0);

% gmm approach
% tame tail; ensure it's below metal band - not necessary
% cap = min(prctile(w, 99.9), thr*0.9);
% w(w > cap) = cap;


% detection approach, for bone
useKDE = true;
try
    [f, xi] = ksdensity(w, 'NumPoints', 2048);
catch
    [cnt, edges] = histcounts(w, 512, 'Normalization','pdf');
    xi = movmean(edges, 2, 'Endpoints','discard');
    f  = movmean(cnt, 5);
end
f = movmean(f, 10);                        % extra smoothing
fr = max(f) - min(f);                      % density dynamic range
if fr <= 0 || ~isfinite(fr), fr = 1; end

% --- Adaptive prominence search for ≥ 2 peaks ---
% Initial guess: 8% of range; relax by 0.7 until >=2 peaks or floor
prom = 0.1 * fr;              
prom_min = 0.01 * fr;          % floor (1% of range)
relax = 0.7;

% Minimal separation between peaks in intensity units (1–3% range)
minDist = 0.1 * (max(xi) - min(xi));

pk = []; pk_idx = [];
while true
    [pk, pk_idx] = findpeaks(f, xi, ...
        'MinPeakProminence', prom, ...
        'MinPeakDistance',   minDist);
    if numel(pk_idx) >= 2 || prom <= prom_min
        break;
    end
    prom = prom * relax;
end

% --- If still <2 peaks, fallback to Otsu ---
if numel(pk_idx) < 2
    xScaled = (w - min(w)) / max(eps, (max(w) - min(w)));
    bone_thr = graythresh(xScaled) * (max(w) - min(w)) + min(w);
    skull_peak_x = NaN;   % for plotting label guard
else
    % Rightmost peak = skull
    skull_peak_x = pk_idx(end);
    prev_peaks = pk_idx(pk_idx < skull_peak_x);

    if ~isempty(prev_peaks)
        left_peak_x = prev_peaks(end);
    else
        % If no earlier peak, expand left until first valley found
        left_peak_x = xi(1);
    end
    
    % Work strictly inside (exclude endpoints)
    in = (xi > left_peak_x) & (xi < skull_peak_x);
    xseg = xi(in); 
    fseg = f(in);
    
    if isempty(xseg)
        % Dead fallback: midpoint
        bone_thr = 0.5*(left_peak_x + skull_peak_x);
    else
        % Robust valley detection in the interval:
        %  absolute minimum value
        fmin = min(fseg);
        %  Tolerance to catch flat valleys/plateaus (fraction of local range)
        frng = max(fseg) - min(fseg);
        tol  = max(1e-12, 1e-3 * max(frng, eps));  % 1e-4 of local range
    
        % Candidate x where f ≈ fmin
        cand = xseg(abs(fseg - fmin) <= tol);
    
        % 4) Remove near-boundary candidates to avoid snapping to edges
        edge_pad = 0.002 * (skull_peak_x - left_peak_x);  % 0.2% of interval
        cand = cand(cand > (left_peak_x + edge_pad) & cand < (skull_peak_x - edge_pad));
    
        if ~isempty(cand)
            % Prefer the RIGHTMOST of the equal-minima candidates
            bone_thr = max(cand);
        else
            % 5) Single argmin with optional quadratic refinement
            [~, j] = min(fseg);
            bone_thr = xseg(j);
    
            % Quadratic interpolation (parabolic fit around j) for sub-bin minimum
            if j > 1 && j < numel(xseg)
                x3 = xseg(j-1:j+1); y3 = fseg(j-1:j+1);
                % Fit y = ax^2 + bx + c; vertex at x* = -b/(2a)
                M = [x3(:).^2, x3(:), ones(3,1)];
                abc = M \ y3(:);
                a = abc(1); b = abc(2);
                if a > 0   % convex
                    xstar = -b/(2*a);
                    if xstar > (left_peak_x + edge_pad) && xstar < (skull_peak_x - edge_pad)
                        bone_thr = xstar;
                    end
                end
            end
        end
    
        % ensure a tiny gap from the skull peak
        min_gap = 0.001 * (max(xi) - min(xi));  % 0.1% of full range
        bone_thr = min(bone_thr, skull_peak_x - min_gap);
    end
end

% against, to speed up, we can use empirical value
if exist('thr','var') && isfinite(thr)
    if ~(bone_thr < thr)
        bone_thr = 0.95 * thr;  % nudge below metal
    end
end

fprintf('Bone thr=%.3f  (prom=%.4g, peaks=%d)\n', bone_thr, prom, numel(pk_idx));

% again disabling plots

assert(bone_thr < thr, 'Bone threshold unexpectedly >= metal threshold.');

BW_skull = (ctVol >= bone_thr) & ~BW;
BW_skull = bwareaopen(BW_skull, 1000);
BW_skull = imclose(BW_skull, strel('sphere', 2));

[X, Y, Z] = meshgrid(1:n, 1:m, 1:p);

% Isosurface - plot disabled



%% 4) Persist into imaging_processing_info_struct and save
% uncomment this for not rewriting
% if ~isfield(imaging_processing_info_struct, 'ct_masks')
%     imaging_processing_info_struct.ct_masks = struct();
% end
ct_masks.metal = struct( ...
    'threshold', thr, ...
    'mask', BW, ...
    'method', 'peak/valley or manual thr', ...
    'notes', 'Cleaned with bwareaopen(500), imclose(sphere,1)' );
ct_masks.skull = struct( ...
    'threshold', bone_thr, ...
    'mask', BW_skull, ...
    'method', 'peak/valley or manual thr', ...
    'notes', 'Cleaned with bwareaopen(1000), imclose(sphere,2)' );



% create struct that stores contact coordinates
% contact_info = struct();
% contact_info stores shank_id, contact_id, contact_type (macro or micro),
% coordinates, and label.
% if use_predefined_location

% this function is for using shank model
 % use_shank_model
    % get shank begin and end
n_trajectory = length(imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory_id);
all_traj_ids = string(imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory_id(:));
n_shank = height(imaging_processing_info_struct.electrode_info);
assert(n_trajectory == n_shank, ...
'Mismatch: %d trajectories vs %d shanks in electrode_info.', n_trajectory, n_shank);
% contact_tbl_path = fullfile(imaging_processing_temp_folder, 'contact_data_table.mat');

% this is only for GT, running locally, for good computational power, just
% reset every time you run it
    % have_contact_tbl = false;
    % if exist(contact_tbl_path, 'file')
    %     Sload = load(contact_tbl_path, 'contact_tbl');
    %     if isfield(Sload, 'contact_tbl') && istable(Sload.contact_tbl)
    %         contact_tbl = Sload.contact_tbl;
    %         have_contact_tbl = true;
    %     end
    % end
    % if ~have_contact_tbl
        % Initialize empty with your exact schema
contact_tbl = table( ...
    'Size',[0 16], ...
    'VariableTypes', {'string','string','uint16','string','double', ...
                      'double','double','double', ...
                      'double','double','double', ...
                      'double','double','double', ...
                      'double','double'}, ...
    'VariableNames', {'Subject','ShankID','ContactID','Model','DeltaMM', ...
                      'X','Y','Z','Ux_lead_dir','Uy_lead_dir','Uz_lead_dir', ...
                       'contact_center_CT_intensity', 'contact_volume_mean_CT_intensity', ...
                       'contact_score',...
                       'radius_mm','length_mm'});

updated_imaging_processing_info_struct_save_path = [temp_folder, ...
        '/imaging_processing_info_struct_after_electrode_localize.mat'];

    
shank_model_all = struct( ...
    'shank_id',{}, 'final_model',{}, 'delta_mm',{}, ...
    'proximal_ras',{}, 'distal_ras',{}, 'bezier_control',{}, ...
    'centerline_pts',{}, 'centerline_s',{}, ...
    'tube_radius_mm',{}, 'tube_auc_line',{}, 'tube_auc_bez',{}, ...
    'BIC_line',{}, 'BIC_bez',{}, 'contact_s_mm' ,{}, ...
    'contact_centers',{}, 'contact_axes',{}, 'shank_n_contact_str',{},...
    'contact_score',{}, 'contact_center_CT_intensity',{}, 'contact_volume_mean_CT_intensity',{},...
    'contact_len_mm',{});

i_traj_list = all_traj_ids;
for i_traj = 1:numel(i_traj_list)
    % i_traj = i_traj_list(i_traj_unprocessed);
    % traj is a string
    
    traj_id = imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory_id{i_traj};
    if opts.Verbose
        fprintf('[%s] Processing shank %d/%d: %s\n', ...
            subj_id, i_traj, numel(i_traj_list), traj_id);
    end
    % check manufacture
    % if ismember(traj_id, string([adtech_research_shank_id{:}]))
    %     current_electrode_manufacture = 'adtech';
    % else
    %     current_electrode_manufacture = electrode_manufacture;
    % end
    current_electrode_manufacture = electrode_manufacture;

    % these are in mm, in RAS coordinates
    start_traj = imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory(i_traj).start;
    end_traj = imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory(i_traj).end;

    % below are known rosa registration erros
    % % !! very important !!
    % if the trajectory in the rosa robot system is wrong, which can
    % happens in the OR, electrode localization will be affected
    % systematic check of correct trajectory is not enabled in this version
    % therefore, it is the user's responsibility to ensure that the end_traj and
    % start_traj are correct

    
    % if strcmp(subj_id, 'BJH046') && strcmp(traj_id, "G'")
    %     end_traj = [18.15, 42.95, 6.41];
    % elseif strcmp(subj_id, 'BJH046') && strcmp(traj_id, "N'")
    %     end_traj = [8.31, 18.80, 39.90];
    % elseif strcmp(subj_id, 'BJH041') && strcmp(traj_id, "A' AMY")
    %     end_traj = [-12.28, -5.75, -24.71];
    % elseif strcmp(subj_id, 'BJH058') && strcmp(traj_id, "O (preCe6 FEF8d prSMA)")
    %     start_traj = [36.00, -10.53, 63.07];
    %     end_traj = [3.53, -9.64, 67.45];
    % elseif strcmp(subj_id, 'BJH055') && strcmp(traj_id, "L'")
    %     end_traj = [-9.38, 23.91, 37.85];
    %     % lead twist (hit hard tissue)
    % elseif strcmp(subj_id, 'BJH040') && strcmp(traj_id, "O PORB OFG")
    %     end_traj = [8.91, 7.98, 11.34];
    % 
    % end

    % viz
    [m,n,p] = size(BW);
    [X, Y, Z] = meshgrid(1:n, 1:m, 1:p);
    % delete video here

    % get the entry point 
    % "proximal to distal" refers to the location
    % of the electrode contacts along the depth of the SEEG electrode
    % as it is inserted into the brain. "Proximal" refers to the
    % electrode contacts closer to the skull and brain surface, while
    % "distal" refers to those deeper within the brain tissue.

    % Sample along line in fine steps (e.g., 0.2 mm)
    % clear cache
    % close all;

    P = sample_line(start_traj, end_traj, 0.2);    % [N x 3]
    ijk1 = (Avox2ras1 \ [P, ones(size(P,1),1)].').'; % 1-based voxel (homog)
    ijk0 = ijk1(:,1:3) - 1;   
    IJ1 = round(ijk0 + 1);            % 1-based indices

    % Clamp to volume bounds
    sz = size(ctVol);
    IJ1 = max(IJ1, 1);
    IJ1(:,1) = min(IJ1(:,1), sz(1));
    IJ1(:,2) = min(IJ1(:,2), sz(2));
    IJ1(:,3) = min(IJ1(:,3), sz(3));

    % Find first skull hit
    linInd = sub2ind(sz, IJ1(:,1), IJ1(:,2), IJ1(:,3));
    hit = find(BW_skull(linInd), 1, 'first');
    if isempty(hit)
        warning('NOT FATAL: No skull intersection found on traj %s; check boneThr.', traj_id);
        entry_ras = start_traj; % use trajectory start point as fallback
    else
        entry_ras = P(hit,:); % proximal (on skull)
    end

% instruction: for correct traj, this parameter does not matter, for
% users using only rosa executed files, this needs to adjust, because there
% are coregistration error between physical space and rosa file space
    roi_rad_mm = 4;
    % 41 A' needs to be special, because two leads are very close
    % if strcmp(subj_id, 'BJH062') && strcmp(traj_id, "E'") ||...
    %      strcmp(subj_id, 'BJH062') && strcmp(traj_id, "G'")||...
    %      strcmp(subj_id, 'BJH056') && strcmp(traj_id, "G") ||...
    %      strcmp(subj_id, 'SLCH024') && strcmp(traj_id, "M") ||...
    %      strcmp(subj_id, 'SLCH020') && strcmp(traj_id, "I") ||...
    %      strcmp(subj_id, 'SLCH017') && strcmp(traj_id, "H SENS SUP")||...
    %      strcmp(subj_id, 'BJH045') && strcmp(traj_id, "N")||...
    %      strcmp(subj_id, 'BJH037') && strcmp(traj_id, "F R ITG FUG EC BENKE")||...
    %      strcmp(subj_id, 'BJH037') && strcmp(traj_id, "B R MTG AHC BENKE")||...
    %      strcmp(subj_id, 'BJH040') && strcmp(traj_id, "Y Y' MFG BIL R CING")
    %     roi_rad_mm = 6;
    % % this is for subject where planing is very far away from
    %     % executed
    % elseif strcmp(subj_id, 'SLCH018') && strcmp(traj_id, "L ORB FRNT CORTEX") ||...
    %      strcmp(subj_id, 'SLCH018') && strcmp(traj_id, "L PST CING") ||...
    %      strcmp(subj_id, 'BJH045') && strcmp(traj_id, "L") ||...
    %      strcmp(subj_id, 'BJH045') && strcmp(traj_id, "M") ||...
    %      strcmp(subj_id, 'BJH045') && strcmp(traj_id, "B") ||...
    %      strcmp(subj_id, 'BJH045') && strcmp(traj_id, "I") ||...
    %      strcmp(subj_id, 'BJH035') && strcmp(traj_id, "E") ||...
    %      strcmp(subj_id, 'BJH040') && strcmp(traj_id, "M PO PERC DLNS ANT")||...
    %      strcmp(subj_id, 'BJH040') && strcmp(traj_id, "O PORB OFG")||...
    %      strcmp(subj_id, 'BJH040') && strcmp(traj_id, "I' ASTG TPC")||...
    %      strcmp(subj_id, 'BJH040') && strcmp(traj_id, "I ITG EC")
    %     roi_rad_mm = 7;
    % elseif strcmp(subj_id, 'BJH045') && strcmp(traj_id, "O") ||...
    %      strcmp(subj_id, 'BJH035') && strcmp(traj_id, "D") 
    %     roi_rad_mm = 8;
    % end

    ROI = cylinder_roi_mask(size(ctVol), Avox2ras0, start_traj, end_traj, roi_rad_mm);
    
    % Quick sanity plot: ROI + planned line (both in RAS) - skipped
   

    % model the shank
    % solution, the entry point (proximal) should be at the skull
    % end point should be around planned trajectory, we will get the
    % end point through connected volume, until reach the end (furthest from the entry point)
    % step_mm = 0.2;
    % connected component might fail because we dont have exact
    % knowledge of the number of components, we can use CC as plan b

    U = end_traj - start_traj;
    L = norm(U);  U = U ./ max(L,eps);
    
    % Tuneables
    
    hi_thr      = thr;                           % your metal peak threshold
    lo_thr      = max(hi_thr - 300, 500);        % lower threshold
    topK         = 6;          % how many seeds to keep (try 4–10)
    min_gap_mm   = 5.0;        % enforce spacing between seeds in mm
    sigma_vox    = 0.6;        % smoothing for peak picking (in voxels)
    dmax_mm = roi_rad_mm / 1.5;
    wI          = 0.7;         % weight: brightness
    wD          = 0.3;         % weight: distance (closer is better)

    % if traj is incorrect, an optional fix is changed lo_thr
    % if  strcmp(subj_id, 'SLCH020') && strcmp(traj_id, "I")
    %     lo_thr = 1000;
    % end

    
    % Masks inside cylinder ROI
    metal_roi_hi = (ctVol >= hi_thr) & ROI;
    metal_roi_lo = (ctVol >= lo_thr) & ROI;
    
    % Local maxima for candidates
    ct_s   = imgaussfilt3(double(ctVol), sigma_vox);
    pk_mask = imregionalmax(ct_s) & metal_roi_lo;
    
    [iPk, jPk, kPk] = ind2sub(size(ctVol), find(pk_mask));
    vals  = ct_s(pk_mask);
    
    % voxel(1-based) -> RAS using 0-based affine
    toRAS = @(i,j,k) ([j-1, i-1, k-1, 1] * Avox2ras0.');  % 1x4
    
    % Distance to the planned line (in RAS, mm)
    cand_ras = zeros(numel(iPk),3);
    rad_mm   = zeros(numel(iPk),1);
    for t = 1:numel(iPk)
        P4 = toRAS(iPk(t), jPk(t), kPk(t));  Pr = P4(1:3);
        cand_ras(t,:) = Pr;
        rad_mm(t) = norm(cross(Pr - start_traj, U));  % perpendicular distance to line
    end
    
    % Discard very far peaks (likely another lead / noise) - optional,
    % we dont know
    % keep = rad_mm <= dmax_mm;
    % iPk = iPk(keep); jPk = jPk(keep); kPk = kPk(keep);
    % vals = vals(keep); cand_ras = cand_ras(keep,:); rad_mm = rad_mm(keep);
   
    
   

    % growing algorithm
    % we have multiple seeds, we will connect the neighbors,
    % to only connect voxels corresponding a lead, we will limit
    % growing direction. In practice, we will continue to grow, and
    % along the process, we will get a tube following the planned
    % trajectory direction that encloses all the connected voxels.

    % fast implementation
    metal_roi = BW & ROI;
    metal_roi_unique = metal_roi;
    % 
  
    % add a assertion here, unique metal roi should be not less than
    % half of the bolb initially detected, for future use - in case
    % modification on growing algo
    assert(( nnz(metal_roi_unique) / nnz(metal_roi_lo) ) > 0.5, ['identified lead is less than half' ...
        ' of the size of metal surroundig planned traj'])

    % assess global directional coherence - function provided - robust
    % implementation ongoing ( parameter learning)
    [ijk_x, ijk_y, ijk_z] = ind2sub(size(metal_roi_unique), find(metal_roi_unique));
    ijk = [ijk_x, ijk_y, ijk_z];  % Nx3 in voxel indices

    opt.CoherencePercent = 50;
    opt.CosineThresh_percentile = 20;
    opt.CosineThresh = 0.98;
    % this changes because the planned trajectory is far away from the
    % roi
    % if strcmp(subj_id, 'SLCH020') && strcmp(traj_id, "I")
    %     opt.CosineThresh = 0.85;
    % end
    % we can also do this
    % keep_mask = cos_vals >= opt.CosineThresh;

    % Transform to RAS coordinates
    P_all = ijk1_to_ras(Avox2ras0, [ijk_y, ijk_x, ijk_z]);
    cosine_stats = voxel_cosine_consistency(P_all, opt.CoherencePercent);
    cos_vals = cosine_stats.avg_cosine_similarity;
    thresh_val = prctile(cos_vals, opt.CosineThresh_percentile);
    % keep_mask = cos_vals >= thresh_val;
    % we do absolute thresholding because depending on the lead, the
    % percentage of voxels corresponding to other lead is different
    keep_mask = cos_vals >= opt.CosineThresh;

    ijk_keep = ijk(keep_mask, :);
    lin_idx_keep = sub2ind(size(metal_roi_unique), ijk_keep(:,1), ijk_keep(:,2), ijk_keep(:,3));
    
    % Build new mask (same size as volume)
    metal_roi_unique_refined = false(size(metal_roi_unique));
    metal_roi_unique_refined(lin_idx_keep) = true;
    metal_roi_unique = metal_roi_unique_refined;
    clear metal_roi_unique_refined;
    
    cmap = turbo;
    cvals = cos_vals;
    


    % remove tiny islands - skip
    

    % now we identify metal volume corresponding to the shank
    % we will define shank entry and end pt, and potentialy the third
    % point, the structure will be added to     
    % shank_models = struct('shank_id',{}, 'proximal_ras',{}, 
    % 'distal_ras',{}, 'dir_unit',{}, 'third_flex_pt',{});
    % we will start with linear fitting as a prior, than do a 3
   % point-intepolation spilne see if it fits the vlume better
    %, approach for quantifying better could be 'baysien factor',
    % the start pt is the one close to the skull
    default_bm = "reference_brainmask.auto.mgz";
    bm_is_default = strcmp(opts.bm_path, default_bm);
    if opts.use_individual_brainmask
        % User requires individual brain mask
        if bm_is_default
            error(['use_individual_brainmask = true, but no bm_path was provided.\n' ...
                   'Provide subject''s brainmask.']);
        end
        if ~isfile(opts.bm_path)
            error('bm_path was provided but does not exist on disk: %s', opts.bm_path);
        end
    else
        % Using reference mask instead
        warning('[%s] use_individual_brainmask = false → reference brain mask is used. Localization may be suboptimal.', subj_id);
    end

    bm_path = opts.bm_path;
    % cfreesurfer fsfast likes char
    bm_path = char(bm_path);

    % Load brainmask, all functionality is included, I disabled plotting
    % for efficiency
    bm = MRIread(bm_path);
    % this electrode bends
    % for bended electrode, we use 'curved', true, this can be added as
    % function input, but with DIXI electrode, the lead is straight!
    % if strcmp(subj_id, 'BJH040') && strcmp(traj_id, "O PORB OFG")
    %     out = fit_shank_line_from_blob(metal_roi_unique, BW_skull, Avox2ras0, start_traj, end_traj, ...
    %         'BrainMaskMGZ', bm, 'BrainMaskVox2RAS', bm.vox2ras0, 'VisualizeFinal', false, 'VisualizeBands', false, ...
    %         'subj_id', subj_id, 'report_dir', temp_folder, 'traj_id', traj_id, 'curved', true);
    % else
    %     out = fit_shank_line_from_blob(metal_roi_unique, BW_skull, Avox2ras0, start_traj, end_traj, ...
    %         'BrainMaskMGZ', bm, 'BrainMaskVox2RAS', bm.vox2ras0, 'VisualizeFinal', false, 'VisualizeBands', false, ...
    %         'subj_id', subj_id, 'report_dir', temp_folder, 'traj_id', traj_id);
    % end

    out = fit_shank_line_from_blob(metal_roi_unique, BW_skull, Avox2ras0, start_traj, end_traj, ...
            'BrainMaskMGZ', bm, 'BrainMaskVox2RAS', bm.vox2ras0, 'VisualizeFinal', false, 'VisualizeBands', false, ...
            'VisualizeTube_rmm', false, 'VisualizeTube_rmm_975', false,...
            'subj_id', subj_id, 'report_dir', temp_folder, 'traj_id', traj_id);

    % when fitting the bolb, we only consider the part within the skull
    % below code is for diagnostic 
   

    if strcmp(current_electrode_manufacture, 'DIXI')
        contact_radius_mm = DIXI_am_diameter/2;
        n_contact_str = imaging_processing_info_struct.electrode_info.n_contact_str{i_traj}; % if exists
        n_contact_num = imaging_processing_info_struct.electrode_info.n_contact(i_traj);     % if exists
        mdl = dixi_model_from_row(n_contact_str, n_contact_num);
    elseif strcmp(current_electrode_manufacture, 'adtech')
        contact_radius_mm = adtech_BF_diameter/2;
        n_contact_str = imaging_processing_info_struct.electrode_info.n_contact_str{i_traj}; % if exists
        n_contact_num = imaging_processing_info_struct.electrode_info.n_contact(i_traj);     % if exists
        mdl = adtech_model_from_row(n_contact_str, n_contact_num);
    elseif strcmp(current_electrode_manufacture, 'PMT')
        contact_radius_mm = PMT_diameter/2;
        n_contact_str = imaging_processing_info_struct.electrode_info.n_contact_str{i_traj}; % if exists
        n_contact_num = imaging_processing_info_struct.electrode_info.n_contact(i_traj);     % if exists
        mdl = pmt_model_from_row(n_contact_str, n_contact_num);
    end

    if strcmp(out.model,'polyline')
        Cpts = out.centerline_pts;
    elseif strcmpi(out.model, 'bezier') && ~isempty(out.bezier_control)
        C_bez = out.bezier_control;
        % Sample densely for arc-length mapping
        [~, Cpts] = bezier_samples(out.proximal_entry_ras, C_bez, out.distal_ras, 1000);
    else
        % Straight line
        tgrid = linspace(0,1,1000).';
        Cpts  = out.proximal_entry_ras + tgrid.*(out.distal_ras - out.proximal_entry_ras);
    end
    seg = vecnorm(diff(Cpts,1,1),2,2);
    s   = [0; cumsum(seg)];            % s=0 at proximal, s(end) at distal
    Ltot = s(end);
    
    % Interp helpers: s -> point, unit tangent
    % pt_at_s = @(ss) [interp1(s, Cpts(:,1), ss, 'linear','extrap'), ...
    %                  interp1(s, Cpts(:,2), ss, 'linear','extrap'), ...
    %                  interp1(s, Cpts(:,3), ss, 'linear','extrap')];
    pt_at_s  = @(ss) point_at_s(out.model, out.proximal_entry_ras, out.distal_ras, ...
        out.bezier_control, out.centerline_pts, ss);


    rad_mm = mdl.diameter_mm/2; 
    len_mm = mdl.contact_len_mm;
    offs   = mdl.offsets_mm(:)';        % from distal tip toward proximal (mm)
    % for MME, we will keep 11 contacts for now
    
    % “s” location of distal tip on the centerline:
    % By construction above, s=0 at proximal (E) and s=Ltot at distal (D).
    s_distal = Ltot;
     % distal tip (deepest point inside brain).
    vx = [imaging_processing_info_struct.ct_volume_registered_to_mri.ysize, ...
              imaging_processing_info_struct.ct_volume_registered_to_mri.xsize, ...
              imaging_processing_info_struct.ct_volume_registered_to_mri.zsize];

    % robust method
    insulating_spacer_length = DIXI_am_spacing - DIXI_am_contact_length;
    score_pose = @(C,u) score_pose_contrast_whole_cylinder_shell(C,u,out.rmm_975,len_mm,Avox2ras0,double(ctVol), vx, insulating_spacer_length);
    % C is the center, u is the direction vector
    ds = min(vx)/2;                                % mm step
    offs_span = max(offs) - min(offs);
    diff_btw_offs_span_and_ltot = max(0, (offs_span - Ltot));
    % s_grid's length should be offs_span + DIXI_am_spacing * 4 / 3
    
    s_grid = ((0 - DIXI_am_spacing * 2 / 3 - len_mm/2 - diff_btw_offs_span_and_ltot):ds:(Ltot + DIXI_am_spacing * 2 / 3)).';
    if diff_btw_offs_span_and_ltot > 0
        assert(abs((max(s_grid) - min(s_grid)) - (offs_span + DIXI_am_spacing * 4 / 3 + len_mm/2)) < ds, 'check s_grid length');
    end
    % we make s_grid larger towards to proximal side, because we the
    % off is center location, we will later subtract len_mm/2
    p = zeros(numel(s_grid),1);
    for i = 1:numel(s_grid)
        Ci = pt_at_s(s_grid(i));
        ui = tan_at_s(s_grid(i), s, Cpts);    
        p(i) = score_pose(Ci, ui);
    end

    win_mm = max(10*len_mm, 20);              % 10x contact or >=20 mm
    win_n  = max(3, round(win_mm/ds));
    base   = movmedian(p, win_n, 'omitnan');
    p_det  = p - base;   

    p_of_s = @(ss) interp1(s_grid, p_det, ss, 'linear', -Inf);
    % interp1(...,'linear', -Inf): linear interpolation; if ss falls outside [min(s_grid), max(s_grid)], return -Inf (your chosen out‑of‑range fill value).

   
% again all intermidiate plotting disabled
    % Intersect with the desired search range, we set it be -2/3
    % spaceing to 1/2 spacing, the spacing is center-to-center distance
    delta_lo = -DIXI_am_spacing * 2 / 3;
    delta_hi = DIXI_am_spacing/2;

    % below adjustment is for cases when the star-traj and end-traj has registration error, 
    % if (strcmp(subj_id, 'SLCH024') && strcmp(traj_id, "E")) ||...
    %    (strcmp(subj_id, 'SLCH024') && strcmp(traj_id, "M")) ||...
    %    (strcmp(subj_id, 'SLCH018') && strcmp(traj_id, "L PST CING"))
    %     delta_lo = -DIXI_am_spacing * 1 / 3;
    %     delta_hi = DIXI_am_spacing * 2 / 3;
    % end
    delta_grid = linspace(delta_lo, delta_hi, 121);

    % delta_grid = linspace(-DIXI_am_spacing, +DIXI_am_spacing, 121);   % widen if needed, in mm
    best_delta = 0; best_obj = -Inf;
    
    for dlt = delta_grid
        % expected centers from distal tip with shift δ
        % Positive δ pushes centers distally (toward larger s).
        % 	With δ=0, the distal-most center (offs=0) sits at s = s_distal − len_mm/2 
        % (i.e., half a contact length shy of the tip
        % s distal is higher
        s_exp = s_distal - (offs(:) + len_mm/2 - dlt);  % column vector
        % if any(s_exp < 0 | s_exp > Ltot), continue; end
        obj = sum(p_of_s(s_exp));                        % sum detrended scores at expected centers
        if obj > best_obj
            best_obj = obj; best_delta = dlt;
        end
    end

% robust method - we show there wont be much difference
    F = griddedInterpolant(s_grid, p_det, 'linear', 'nearest');     % fast p_det(s)
    p_of_s_safe = @(ss) deal_vals(ss, F, min(s_grid), max(s_grid));  
    % diagnostic
    % fprintf('F covers [%.2f, %.2f]\n', min(s_grid), max(s_grid));
    % fprintf('delta_lo = %.2f, Ltot = %.2f, delta_hi = %.2f\n', delta_lo, Ltot, delta_hi);

    
    % Expected contact centers given a shift dlt
    % if dlt = 0, the smallest exp_center will be at len_mm/2, and
    % lower dlt (for example -1), means the exp center is shifted
    % towards proximal side
    exp_centers = @(dlt) s_distal - (offs(:) + len_mm/2 - dlt);     % column vector
    % s_test = exp_centers(0);  % delta = 0
    % fprintf('exp_centers(0): min = %.2f, max = %.2f\n', min(s_test), max(s_test));
    % s_test = exp_centers(delta_grid(1));  % delta = 0
    % fprintf('exp_centers(delta_grid(1)): min = %.2f, max = %.2f\n', min(s_test), max(s_test));
    
    % Robust, trimmed-sum objective (tolerates missing/weak contacts)
    trim_q = 0.80;                                                  % keep top 80% contacts
    robj = @(dlt) robust_sum(p_of_s_safe(exp_centers(dlt)), trim_q);
    
    obj_grid_orig = arrayfun(@(dlt) sum(p_of_s_safe(exp_centers(dlt))), delta_grid);
    obj_grid_rob  = arrayfun(@(dlt)             robj(dlt)         , delta_grid);
    
    [~,ix_orig] = max(obj_grid_orig);
    best_delta_orig = delta_grid(ix_orig);
    
    [~,ix_rob]  = max(obj_grid_rob);
    best_delta_rob0 = delta_grid(ix_rob);
    
    % Narrow continuous refine for the robust objective
    lb = max(delta_grid(1), best_delta_rob0 - 0.25);
    ub = min(delta_grid(end), best_delta_rob0 + 0.25);
    best_delta_rob = fminbnd(@(d)-robj(d), lb, ub);                 % maximize
    
    % diagnostic plot 
    % Expected centerlines from both deltas
    s_exp_orig = exp_centers(best_delta_orig);
    s_exp_rob  = exp_centers(best_delta_rob);
    
   

    delta_use = best_delta_orig;      % or: best_delta_rob
    % delta_use = best_delta_rob;

    % Expected centers from the chosen global delta (comb)
    s_seed = exp_centers(delta_use);                 % N×1
    % no need to clip it, if electrode is outside of skull, let it be
    % s_seed = min(max(s_seed,0), Ltot + DIXI_am_spacing / 2);               % clip to domain
    
    % --- PER-CONTACT MICRO REFINE (around each seed)
    ds     = s_grid(2)-s_grid(1);
    ax_win = DIXI_am_spacing / 20; 

    % adjustment needed when start traj and end traj are not correct
    % if strcmp(subj_id, 'BJH040') && strcmp(traj_id, "O PORB OFG")
    %     ax_win = DIXI_am_spacing / 5;
    % end
    
    s_refined = nan(size(s_seed));
    for k = 1:numel(s_seed)
        s0 = s_seed(k);
        lb = max(0,    s0 - ax_win);
        ub = min(Ltot, s0 + ax_win);
        if ub > lb
            s_refined(k) = fminbnd(@(ss) -p_of_s_safe(ss), lb, ub);  % maximize p_det
        else
            s_refined(k) = s0;   % fallback
        end
    end
    % nexttile(2); % second panel in the tiledlayout
    s_refined_orig = s_exp_orig;
    % s_refined_orig = s_exp_rob;
    
    if strcmp(mdl.type, 'MME')
        contact_s = s_refined(mdl.active_idx);
    else
        contact_s = s_refined(:);         
    end
    % Nx1
    % stem_safe(contact_s, p_of_s_safe(contact_s), 'Marker','o','Color',[0.05 0.90 0], 'DisplayName','refined (orig)');

    cellsC = arrayfun(@(ss) pt_at_s(ss), contact_s, 'UniformOutput', false);
    contact_centers = cat(1, cellsC{:});     % or: vertcat(cellsC{:})
    
    % Tangents (N×3)
    cellsU = arrayfun(@(ss) tan_at_s(ss, s, Cpts), contact_s, 'UniformOutput', false);
    contact_tangents = cat(1, cellsU{:});
    
    % along-track score (N×1)
    contact_score = p_of_s_safe(contact_s);
    
    % CT intensities (nearest voxel + small orthogonal disc average),
    % in sample_ct_nn we must do ijk index
    contact_HU_nn   = sample_ct_nn(double(ctVol), Avox2ras0, contact_centers);           % N×1
    contact_HU_disc = arrayfun(@(k) meanHU_disc(double(ctVol), Avox2ras0, ...
                                   contact_centers(k,:), contact_tangents(k,:), ...
                                   0.6*(mdl.diameter_mm/2), 0.6*len_mm), (1:numel(contact_s))');


    
    % reorganize output to struct
    % put electrode coordinates in to contact_data_table
    fn = fieldnames(shank_model_all);
    shank_model = cell2struct(repmat({[]}, 1, numel(fn)), fn, 2);
    
    % Required fields
    shank_model.shank_id      = traj_id;
    shank_model.proximal_ras  = out.proximal_entry_ras(:).';
    shank_model.distal_ras    = out.distal_ras(:).';
    shank_model.centerline_pts= Cpts;
    shank_model.centerline_s  = s;
    shank_model.tube_radius_mm= rad_mm;
    shank_model.delta_mm      = best_delta_orig;
    

    shank_model.final_model          = out.model;
    shank_model.bezier_control = out.bezier_control;  % 1x3

    
    % Scores (use safe getters; fall back to NaN if not present)
    shank_model.tube_auc_line = getf(out,'tube_auc_line', NaN);
    shank_model.tube_auc_bez  = getf(out,'tube_auc_bezier',  NaN);   % ensure your upstream uses this exact name
    shank_model.BIC_line      = getf(out,'bic_line',      NaN);   % match case: 'BIC_line'
    shank_model.BIC_bez       = getf(out,'bic_bezier',       NaN);
    
    
    % Contacts
    shank_model.contact_s_mm     = contact_s(:);
    shank_model.contact_centers   = contact_centers;
    shank_model.contact_axes      = contact_tangents;
    shank_model.contact_score    = contact_score(:);             % p_det at s
    shank_model.contact_center_CT_intensity    = contact_HU_nn(:);             % nearest-voxel HU
    shank_model.contact_volume_mean_CT_intensity  = contact_HU_disc(:);  
    shank_model.contact_len_mm   = len_mm;
    shank_model.shank_n_contact_str = n_contact_str;
    
    % 
    shank_model_all(end+1) = shank_model;

    % for k = 1:size(contact_centers,1)
    %     Ck = contact_centers(k,:); uk = contact_tangents(k,:);
    %     contacts_cell(end+1,:) = {subj_id, traj_id, k, shank_model.final_model, best_delta, ...
    %                               Ck(1),Ck(2),Ck(3), uk(1),uk(2),uk(3), ...
    %                               contact_intensity(k), rad_mm, len_mm};
    % end

    nK = size(contact_centers, 1);

    Tk = table( ...
        repmat(string(subj_id), nK, 1), ...
        repmat(string(traj_id), nK, 1), ...
        uint16((1:nK).'), ...
        repmat(string(shank_model.final_model), nK, 1), ...
        repmat(best_delta_orig, nK, 1), ...
        contact_centers(:,1), contact_centers(:,2), contact_centers(:,3), ...
        contact_tangents(:,1), contact_tangents(:,2), contact_tangents(:,3), ...
        contact_HU_nn, contact_HU_disc, contact_score,...
        repmat(rad_mm, nK, 1), ...
        repmat(len_mm, nK, 1), ...
        'VariableNames', contact_tbl.Properties.VariableNames);
    % MeanCT – Mean CT intensity inside the cylindrical ROI for this contact (proxy for metal overlap).
    % RadMM – Electrode contact radius in mm.
    % LenMM – Electrode contact length in mm.
    % Ux, Uy, Uz – Components of the unit tangent vector along the electrode at this contact’s center (RAS coordinates).
    % Defines the cylinder’s axis direction for intensity sampling and visualization.
    % ContactID – Index of the contact along this shank (1 = most distal, increasing proximally).
    contact_tbl = [contact_tbl; Tk];

    % save as we go, as this process might take about 10 minutes per
    % lead
    table_filename = fullfile(temp_folder, [subj_id,  ...
        '_contact_data_table.mat']);
    save(table_filename, 'contact_tbl', '-v7.3');
    contact_data_struct = table2struct(contact_tbl);
    table_filename = fullfile(temp_folder, ...
        [subj_id, '_contact_data_table_struct.mat']);
    save(table_filename, 'contact_data_struct', '-v7.3');
    imaging_processing_info_struct.shank_model_all = shank_model_all;
    save(fullfile(temp_folder,  [subj_id, '_imaging_processing_info_struct_update.mat']), 'imaging_processing_info_struct', '-v7.3');

end

        

% export the realigned electrode location and check in freeview
unique_shankIDs = unique(string(contact_tbl.ShankID));
for i_shankID = 1:length(unique_shankIDs)
    shankID = unique_shankIDs{i_shankID};
    filtered_table = contact_tbl(strcmp(contact_tbl.ShankID, shankID), :);

    filtered_table = sortrows(filtered_table, 'ContactID', 'ascend');

    % Assemble XYZ (RAS, mm)
    XYZ = [filtered_table.X, filtered_table.Y, filtered_table.Z];
    file_name = sprintf('%s.dat', shankID);
    % fid = fopen([automated_aligned_to_brightest_voxel_macro_contacts_folder, '/', file_name], 'w');
    fid = fopen(fullfile(automated_aligned_to_brightest_voxel_macro_contacts_folder, file_name), 'w');
    assert(fid>0, 'Cannot open %s for writing.', file_name);
    fprintf(fid, '\n');
    fprintf(fid, '%.6f %.6f %.6f\n', XYZ.');

    % for i_contact = 1:height(filtered_table)
    %     fprintf(fid, [num2str([X(i_contact), Y(i_contact), Z(i_contact)]) '\n']);
    % end
    fprintf(fid, 'info\n');
    fprintf(fid, 'numpoints %d\n', height(filtered_table));
    fprintf(fid, 'useRealRAS 1\n');
    fclose(fid);
end

% for diagnostic
for i = 1:numel(shank_model_all)
    if isempty(shank_model_all(i).centerline_pts), continue; end
    C = shank_model_all(i).centerline_pts;
    file_name = sprintf('%s.dat', shank_model_all(i).shank_id);
    % fid = fopen([automated_shank_center_line_folder, '/', file_name], 'w');
    fid = fopen(fullfile(automated_shank_center_line_folder, file_name), 'w');
    fprintf(fid, '%.6f %.6f %.6f\n', C.'); 
    fprintf(fid, 'info\nnumpoints %d\nuseRealRAS 1\n', size(C,1)); fclose(fid);
end



%% Micro wires localization, contact labeling and segmentation are not implemented in this function
%% but will be included in the full package.

