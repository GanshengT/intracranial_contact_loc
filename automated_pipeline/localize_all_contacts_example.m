%% Description
% This script localize the micro contact, adjust manually localized macro
% contact
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
% transformation of the electrode to get brain parcellation label. Also, we
% usually do not reslice ct for electrode localization.

% contact: g.tan@wustl.edu
clear;
close all;
%% definition
% we load all the info, add all the dependent function
% for example, addpath(genpath(..... script/imaging_processing_scripts']));
data_path = '...';
% subj info
subj_id = 'BJH079';
electrode_manufacture = 'DIXI';
DIXI_am_diameter = 0.8;
DIXI_am_contact_length = 2; % mm
DIXI_am_spacing = 3.5;
DIXI_cm_larger_spacing = 9;
DIXI_bm_larger_spacing = 13;
use_predefined_location = false;
use_shank_model = true;
% research shank name
% this is subject-specific, be sure that you change that. This should match
% manually identified labels
research_shank_id = {'HL_MTG_AHC_BEFR', 'HR_MTG_AHC_BEFR'};
% no research electrode
research_shank_id = {};
% manufacture spect, might be useful
elec_model_dixi_15cm = [0 3.5 7 10.5 14 27 30.5 34 37.5 41 54 57.5 61 64.5 68];
elec_model_dixi_18cm = [0 3.5 7 10.5 14 17.5 30.5 34 37.5 41 44.5 48 61 64.5 68 71.5 75 78.5];
elec_model_dixi_am = 0:3.5:59.5;
%elec_model_adtech = 0:5:35;
elec_model_dixi_mme   = [0 4 8 32 36 40];

% research_shank_id = {'A''^MTG-sts-Amy', 'B^MTG-aHc'};
% make the string array into char array
% research_shank_id = cellfun(@char, research_shank_id, 'UniformOutput', false);
% based on ad-tech catalog
% the macro contact is 1.57mm long, about 1-3 pixels
% with 1.28mm diameter, 1-2 pixels
% we will fit gaussian on a subvolume in CT centered at each manually
% localizated macro electrode
sub_vol_range_mm = 1.57; 
% mm unit, we want the subvol engulge the contact
interp_ratio = 10; 
% this parameter is no longer used
% to adjust macro contact location to the brightest voxel (most likely
% contact location), we will break a voxel into 10*10*10 small voxels.
% the current method to get objective macro electrode location:
% gradually decrease thresholding intensity until we found a connected
% volume
% continuous increase the volume if the connected volume does not meet the
% minimum_volume_for_contact
intensity_search_resolution = 1;
macro_contact_diameter = 1.28;
macro_contact_length = 1.57;
macro_contact_volume_mm_cube = macro_contact_diameter^2 * pi * macro_contact_length;
micro_contact_from_first_macro = 4; 
minimum_voxel_for_micro_contact = 2; % this is empirical, adjust if needed
% we assume micro contact is 3mm from the first macro
isovalue_to_get_surface = 0.5;
% 0.5 is chosen for getting surface for binary volume
% we will use nearest neighbor for label assignment, we first define the
% volume for neighbor to be considered
neighbor_distance_requirement = sqrt(macro_contact_diameter^2 + ...
    macro_contact_length^2); % mm
shank_start_from_last_macro_contact = 10; %mm, for plotting shank purpose
shank_tip_diameter = 0.2;
shank_tip_length = 1;
% the tip is a cylinder of 0.2mm diameter of 1mm length
micro_contact_id = 99;
imaging_processing_temp_folder = [data_path, '/', subj_id, '/imaging_process'];
imaging_processing_info_struct_save_path = [imaging_processing_temp_folder, ...
    '/imaging_processing_info_struct.mat'];
imaging_processing_info_struct = load(imaging_processing_info_struct_save_path);

% if you do not have manually identified macro contact ready yet - option
% one
% uncomment below for manual edition in freeview
% !! important, from proximal to distal, id number goes from 1 to 16 !!
% OR_sheet_path = [data_path, '/' subj_id, '/notes/SEEG-OR-', subj_id, '.xlsx'];
% sheetname = 'Surgeon Date';
% col_name_electrode_number = 'Est_Electrode';
% electrode_info = ...
% extract_electrode_info_fromOR_sheet(OR_sheet_path, sheetname, col_name_electrode_number);


% update note
% we will identify shank from CT volume and trajectory derived from ROSA
% so at the end, we will have estimated anchor points for the shank, the
% shank is model by intepolation of 2 or more anchor pts
% During this process, we will estimate the location for each macro
% contacts, specifically, we do objective function optimization. Based on
% electrode distance, electrode number, we find the points along modeled
% shank. 

% Objective function, modeled shank and electrode vs bright voxels.
% How to do optimization, estimate the first electrode location and
% distance


% load manually identified macro contact
manually_identified_macro_contacts_path = [imaging_processing_temp_folder, ...
    '/manually_identified_macro_contacts/'];
% save macro contacts' location that are aligned with the brightest voxel
aligned_to_brightest_voxel_macro_contacts_folder = [imaging_processing_temp_folder, ...
    '/aligned_to_brightest_voxel_macro_contacts'];
if ~exist(aligned_to_brightest_voxel_macro_contacts_folder, 'dir')
    mkdir(aligned_to_brightest_voxel_macro_contacts_folder);
end

automated_aligned_to_brightest_voxel_macro_contacts_folder = [imaging_processing_temp_folder, ...
    '/automated_aligned_to_brightest_voxel_macro_contacts'];
if ~exist(automated_aligned_to_brightest_voxel_macro_contacts_folder, 'dir')
    mkdir(automated_aligned_to_brightest_voxel_macro_contacts_folder);
end

automated_shank_center_line_folder = [imaging_processing_temp_folder, ...
    '/automated_shank_center_lines'];
if ~exist(automated_shank_center_line_folder, 'dir')
    mkdir(automated_shank_center_line_folder);
end

% save contacts' location visualization
report_contacts_related_surf_folder = [imaging_processing_temp_folder, ...
    '/report_contacts_related_surf_folder'];
if ~exist(report_contacts_related_surf_folder, 'dir')
    mkdir(report_contacts_related_surf_folder);
end

% region of interest
% we will focus on aparc - Desikan/Killiany parcellation
region_of_interest = {'Left-Cerebral-White-Matter', 'Left-Thalamus', ...
    'Left-Caudate', 'Left-Putamen', ...
    'Left-Pallidum', 'Left-Hippocampus', 'Left-Amygdala', 'Left-Accumbens-area', ...
    'Left-VentralDC', 'ctx-lh-bankssts', 'ctx-lh-caudalanteriorcingulate', ...
    'ctx-lh-caudalmiddlefrontal', 'ctx-lh-corpuscallosum', 'ctx-lh-cuneus', ...
    'ctx-lh-entorhinal', 'ctx-lh-fusiform', 'ctx-lh-inferiorparietal', ...
    'ctx-lh-inferiortemporal', 'ctx-lh-isthmuscingulate', ...
    'ctx-lh-lateraloccipital', 'ctx-lh-lateralorbitofrontal', ...
    'ctx-lh-lingual', 'ctx-lh-medialorbitofrontal', 'ctx-lh-middletemporal', ...
    'ctx-lh-parahippocampal', 'ctx-lh-paracentral', 'ctx-lh-parsopercularis', ...
    'ctx-lh-parsorbitalis', 'ctx-lh-parstriangularis', 'ctx-lh-pericalcarine', ...
    'ctx-lh-postcentral', 'ctx-lh-posteriorcingulate', 'ctx-lh-precentral', ...
    'ctx-lh-precuneus', 'ctx-lh-rostralanteriorcingulate', ...
    'ctx-lh-rostralmiddlefrontal', 'ctx-lh-superiorfrontal', ...
    'ctx-lh-superiorparietal', 'ctx-lh-superiortemporal', ...
    'ctx-lh-supramarginal', 'ctx-lh-frontalpole', 'ctx-lh-temporalpole', ...
    'ctx-lh-transversetemporal', 'ctx-lh-insula'};
region_of_interest_hippo_subsegment = {...
    'presubiculum-head', 'subiculum-head', 'CA1-head', 'CA3-head', ...
    'CA4-head', 'GC-ML-DG-head', 'molecular_layer_HP-head', 'Lateral-nucleus', ...
    'Basal-nucleus', 'Central-nucleus', 'Medial-nucleus', 'Cortical-nucleus'...
    'Accessory-Basal-nucleus'
    };
    
    % 'ctx_lh_G_and_S_frontomargin', 'ctx_lh_G_and_S_transv_frontopol', ...
    % 'ctx_lh_G_front_inf-Opercular', 'ctx_lh_G_front_inf-Orbital', ...
    % 'ctx_lh_G_front_inf-Triangul', 'ctx_lh_G_front_middle', 'ctx_lh_G_front_sup', ...
    % 'ctx_lh_G_oc-temp_lat-fusifor', 'ctx_lh_G_oc-temp_med-Lingual', ...
    % 'ctx_lh_G_oc-temp_med-Parahip', 'ctx_lh_G_orbital', 'ctx_lh_G_temp_sup-G_T_transv',...
    % 'ctx_lh_G_temp_sup-Lateral', 'ctx_lh_G_temp_sup-Plan_polar', ...
    % 'ctx_lh_G_temp_sup-Plan_tempo', 'ctx_lh_G_temporal_inf', 'ctx_lh_G_temporal_middle',...
    % 'ctx_lh_Pole_temporal', 'ctx_lh_S_front_inf', 'ctx_lh_S_front_middle',...
    % 'ctx_lh_S_front_sup', 'ctx_lh_S_temporal_inf', 'ctx_lh_S_temporal_sup',...
    % 'ctx_lh_S_temporal_transverse',...
    
region_of_interest_to_add = {};
for i_region_of_interest = 1:length(region_of_interest)
    current_region = region_of_interest{i_region_of_interest};
    if contains(current_region, 'Left')
        new_region = strrep(current_region, 'Left', 'Right');
        region_of_interest_to_add = [region_of_interest_to_add, new_region];
    end
    
    if contains(current_region, 'lh')
        new_region = strrep(current_region, 'lh', 'rh');
        region_of_interest_to_add = [region_of_interest_to_add, new_region];
    end
end

region_of_interest = [region_of_interest, region_of_interest_to_add];
region_of_interest_to_add = {};
for i_region_of_interest = 1:length(region_of_interest_hippo_subsegment)
    current_region = region_of_interest_hippo_subsegment{i_region_of_interest};
    if contains(current_region, 'head')
        new_region = strrep(current_region, 'head', 'body');
        region_of_interest_to_add = [region_of_interest_to_add, new_region];
    end
end
region_of_interest = [region_of_interest, region_of_interest_to_add];

% get CT intensity threshold, for bone and electrode volume
ctVol = imaging_processing_info_struct.ct_volume_registered_to_mri.vol;      % 3D numeric
Avox2ras1 = imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras1; % 4x4 affine (voxel->RAS)
Avox2ras0 = imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras0; % 4x4 affine (voxel->RAS)

v = double(ctVol(:));
v = v(isfinite(v));
v = v(v > 0);                      % drop air if desired
[nPDF.x, nPDF.f] = deal([]);
try
    [f, xi] = ksdensity(v, 'NumPoints', 1024);  % smooth density
    nPDF.x = xi; nPDF.f = f;
catch
    % fallback to histogram-based density if KDE unavailable
    [cnt, edges] = histcounts(v, 256, 'Normalization','pdf');
    xi = movmean(edges,2,'Endpoints','discard');
    f  = movmean(cnt,3); % smooth a bit
    nPDF.x = xi; nPDF.f = f;
end

% Find peaks (prominence guards spurious ones)
[pk, locs] = findpeaks(nPDF.f, nPDF.x, 'MinPeakProminence', max(nPDF.f)*0.02);

thr = prctile(v,99);  % default fallback
if numel(locs) >= 2
    % take two highest-x peaks (tissue + metal)
    [~,ix] = sort(locs,'ascend');
    lo = locs(ix(end-1)); hi = locs(ix(end));
    % optionally enforce metal peak > 3000 (adjust to your units)
    if hi > 2900  % your metal guard
        in = nPDF.x >= lo & nPDF.x <= hi;
        xseg = nPDF.x(in);
        fseg = nPDF.f(in);

        % find minimal density value in the interval
        fmin = min(fseg);
        % tolerance to catch flat minima/plateaus
        tol  = max(1e-12, 1e-4 * range(fseg));
        cand = xseg(abs(fseg - fmin) <= tol);

        % take the largest threshold among minima
        if ~isempty(cand)
            thr = max(cand);
        else
            % numeric fallback: just take the single argmin
            [~,j] = min(fseg);
            thr = xseg(j);
        end
    end
end

assert((thr <= 4000) & (thr >= 2000), ...
    'identified thres is not between 2000 - 4000, rule of sum');

% viz to confirm
figure('Color','w'); hold on;
% Histogram as PDF
nbins = 256;
histogram(v, nbins, 'Normalization','pdf', 'EdgeColor','none', 'FaceAlpha',0.35);
% KDE for smooth view
try
    [f, xi] = ksdensity(v, 'NumPoints', 512);
    plot(xi, f, 'LineWidth', 1.5);
end
xline(thr, ':', 'identified thres', 'LabelVerticalAlignment','middle');
title('distribution of ct intensity and threshold for isosurfacing electrodes')
% fill your report path
% report_dir = ['...../report/', subj_id, '/imaging_process'];
if ~exist(report_dir, 'dir')
    mkdir(report_dir);
end

% Generate filename with subject ID and descriptive name
filename = [subj_id, '_CT_intensity_distribution_electrode_threshold'];
full_path = fullfile(report_dir, [filename, '.pdf']);
saveas(gcf, full_path, 'pdf');


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

figure('Color','w'); hold on;
h = patch(isosurface(X, Y, Z, BW, 0.5));             % use meshgrid coords
isonormals(X, Y, Z, double(ctVol), h);               % normals from intensity
set(h, 'EdgeColor','none', 'FaceColor',[0.75 0.78 0.85], 'FaceAlpha',0.25);

% decimate to speed up rendering
% h = reducepatch(h, 0.5);  % keep 50% of faces
% Transform vertices from voxel (1-based) -> 0-based -> RAS (mm)
Vvox1 = h.Vertices;           % [i j k] in 1-based voxel indices
% V0    = Vvox1 - 1;            % 0-based for vox2ras1
% Vras  = ( [V0, ones(size(V0,1),1)] * Avox2ras.' );
Vras  = ([Vvox1, ones(size(Vvox1,1),1)] * Avox2ras1.');
h.Vertices = Vras(:,1:3);


axis equal vis3d; camlight headlight; lighting gouraud;
grid on; box on;
xlabel('X_R_A_S (mm)'); ylabel('Y_R_A_S (mm)'); zlabel('Z_R_A_S (mm)');
title(sprintf('CT ≥ %.0f (electrode isosurface in RAS)', thr));
filename = [subj_id, '_identified_electrode_isosurface'];
full_path = fullfile(report_dir, [filename, '.pdf']);
saveas(gcf, full_path, 'pdf');

% get skull mask
ctNoMetal = double(ctVol);
ctNoMetal(ctVol >= thr) = NaN; 
w = ctNoMetal(:);
w = w(isfinite(w) & w > 0);

% gmm approach
% tame tail; ensure it's below metal band - not necessary
% cap = min(prctile(w, 99.9), thr*0.9);
% w(w > cap) = cap;

% bone_thr = [];
% try
%     gm = fitgmdist(w, 2, 'RegularizationValue', 1e-6, 'Replicates', 3, ...
%         'Options', statset('MaxIter', 300, 'Display', 'off'));
%     % sort components by mean
%     [mus, ord] = sort(gm.mu(:), 'ascend');
%     ps = gm.ComponentProportion(:);  ps = ps(ord);
%     % handle Sigma shape
%     S = gm.Sigma; if ndims(S)==3, S = squeeze(S); end
%     if numel(S)==1, S = [S S]; end
%     s1 = sqrt(S(ord(1))); s2 = sqrt(S(ord(2)));
% 
%     % intersection of two normals: solve p1*N(mu1,s1) = p2*N(mu2,s2)
%     mu1 = mus(1); mu2 = mus(2);
%     p1  = ps(1);  p2  = ps(2);
% 
%     a = 1/(2*s2^2) - 1/(2*s1^2);
%     b = mu1/(s1^2) - mu2/(s2^2);
%     c = (mu2^2)/(2*s2^2) - (mu1^2)/(2*s1^2) - log((p2*s1)/(p1*s2));
%     xs = roots([a b c]);                             % up to two reals
%     xs = xs(isreal(xs));
%     % pick the root between the two means (or nearest to midpoint)
%     if ~isempty(xs)
%         inBetween = xs(xs >= mu1 & xs <= mu2);
%         if ~isempty(inBetween)
%             bone_thr = inBetween(1);
%         else
%             [~,kmin] = min(abs(xs - (mu1+mu2)/2));
%             bone_thr = xs(kmin);
%         end
%     else
%         bone_thr = (mu1 + mu2) / 2;
%     end
% 
% catch
%     % fallback if GMM fails: Otsu on nonmetal intensities
%     xScaled = (w - min(w)) / max(eps, (max(w) - min(w)));
%     bone_thr = graythresh(xScaled) * (max(w) - min(w)) + min(w);
% end

% peak detection approach
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
        % 1) Find absolute minimum value
        fmin = min(fseg);
    
        % 2) Tolerance to catch flat valleys/plateaus (fraction of local range)
        frng = max(fseg) - min(fseg);
        tol  = max(1e-12, 1e-3 * max(frng, eps));  % 1e-4 of local range
    
        % 3) Candidate x where f ≈ fmin
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
    
        % 6) Ensure a tiny gap from the skull peak
        min_gap = 0.001 * (max(xi) - min(xi));  % 0.1% of full range
        bone_thr = min(bone_thr, skull_peak_x - min_gap);
    end
end

% --- Sanity guards relative to metal threshold ---
if exist('thr','var') && isfinite(thr)
    if ~(bone_thr < thr)
        bone_thr = 0.95 * thr;  % nudge below metal
    end
end

fprintf('Bone thr=%.3f  (prom=%.4g, peaks=%d)\n', bone_thr, prom, numel(pk_idx));

% Diagnostics plot ----------
figure('Color','w'); hold on;
histogram(w, 256, 'Normalization','pdf', 'EdgeColor','none', 'FaceAlpha',0.25);
plot(xi, f, 'LineWidth', 1.5);

% mark peaks/valleys if present
if ~isempty(pk_idx)
    plot(pk_idx, pk, 'v', 'MarkerFaceColor',[0.2 0.2 0.2], 'MarkerEdgeColor','none');
end
% [~, v_idx_all] = findpeaks(-f);
% if ~isempty(v_idx_all)
%     plot(xi(v_idx_all), f(v_idx_all), '^', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor','none');
% end

if exist('skull_peak_x','var') && isfinite(skull_peak_x)
    xline(skull_peak_x, '-', 'skull peak', 'LabelVerticalAlignment','bottom');
end
xline(bone_thr, '-', 'bone thr', 'LabelVerticalAlignment','bottom', 'LineWidth',1.2);
if exist('thr','var'), xline(thr, ':', 'metal thr'); end

xlabel('CT intensity'); ylabel('Density');
ttl = sprintf('Peak→Valley bone threshold (MinProm ~ %.1f%% of range)', 100*prom/fr);
title(ttl);
grid on; box on;
filename = [subj_id, '_identified_skull_based_on_intensity'];
full_path = fullfile(report_dir, [filename, '.pdf']);
saveas(gcf, full_path, 'pdf');


assert(bone_thr < thr, 'Bone threshold unexpectedly >= metal threshold.');

BW_skull = (ctVol >= bone_thr) & ~BW;
BW_skull = bwareaopen(BW_skull, 1000);
BW_skull = imclose(BW_skull, strel('sphere', 2));

% gmm fitting
% figure('Color','w'); hold on;
% histogram(w, 256, 'Normalization','pdf', 'EdgeColor','none', 'FaceAlpha',0.35);
% try
%     [f, xi] = ksdensity(w, 'NumPoints', 512); plot(xi, f, 'LineWidth', 1.5); end
% xline(bone_thr, '-', 'bone thr', 'LabelVerticalAlignment','bottom');
% xline(thr,      ':', 'metal thr', 'LabelVerticalAlignment','bottom');
% xplot = linspace(prctile(w,0.1), prctile(w,99.9), 512);
% 
% % Robustly extract/sort component params
% [mus, ord] = sort(gm.mu(:), 'ascend');
% ps = gm.ComponentProportion(:); ps = ps(ord);
% S  = gm.Sigma; if ndims(S)==3, S = squeeze(S); end
% if numel(S)==1, S = [S S]; end
% s1 = sqrt(S(ord(1))); s2 = sqrt(S(ord(2)));
% 
% % Component PDFs and mixture
% g1 = ps(1) * normpdf(xplot, mus(1), s1);
% g2 = ps(2) * normpdf(xplot, mus(2), s2);
% mix = g1 + g2;
% 
% plot(xplot, g1, '--', 'LineWidth', 1.2);
% plot(xplot, g2, '--', 'LineWidth', 1.2);
% plot(xplot, mix, '-',  'LineWidth', 1.5);
% title('CT intensity (non-metal) with bone threshold'); grid on; box on; xlabel('CT intensity'); ylabel('pdf');
% 
% % figure for bone
% [m,n,p] = size(BW_skull);
[X, Y, Z] = meshgrid(1:n, 1:m, 1:p);

% Isosurface
figure('Color','w'); hold on;
hb = patch(isosurface(X, Y, Z, BW_skull, 0.5));
isonormals(X, Y, Z, double(ctVol), hb);
set(hb, 'EdgeColor','none', 'FaceColor',[0.85 0.80 0.70], 'FaceAlpha',0.35); % bone-ish

% Transform vertices to RAS (1-based -> 0-based -> RAS)
Vvox1 = hb.Vertices;            % [x y z] in 1-based voxel indices
% V0    = Vvox1 - 1;              % 0-based for vox2ras1
Vras  = ([Vvox1, ones(size(Vvox1,1),1)] * Avox2ras1.');
hb.Vertices = Vras(:,1:3);

axis equal vis3d; camlight headlight; lighting gouraud;
grid on; box on;
xlabel('X_RAS (mm)'); ylabel('Y_RAS (mm)'); zlabel('Z_RAS (mm)');
title(sprintf('Skull isosurface (CT %.0f, threshold, metal removed)', bone_thr));
filename = [subj_id, '_identified_skull_isosurface'];
full_path = fullfile(report_dir, [filename, '.pdf']);
saveas(gcf, full_path, 'pdf');


%% 4) Persist into imaging_processing_info_struct and save
% uncomment this for not rewriting
% if ~isfield(imaging_processing_info_struct, 'ct_masks')
%     imaging_processing_info_struct.ct_masks = struct();
% end
imaging_processing_info_struct.ct_masks.metal = struct( ...
    'threshold', thr, ...
    'mask', BW, ...
    'method', 'peak/valley or manual thr', ...
    'notes', 'Cleaned with bwareaopen(500), imclose(sphere,1)' );

imaging_processing_info_struct.ct_masks.skull = struct( ...
    'threshold', bone_thr, ...
    'mask', BW_skull, ...
    'method', 'peak/valley or manual thr', ...
    'notes', 'Cleaned with bwareaopen(1000), imclose(sphere,2)' );



% create struct that stores contact coordinates
% contact_info = struct();
% contact_info stores shank_id, contact_id, contact_type (macro or micro),
% coordinates, and label.
if use_predefined_location
    datFiles = dir(fullfile(manually_identified_macro_contacts_path, '*.dat'));
    % Initialize a cell array to store the contact data
    contact_data = cell(0, 5); % Initialize an empty cell array with 5 columns
    
    % Loop through each .dat file
    for i_dat_file = 1:numel(datFiles)
        skipped_line = 0;
        [~, shankID, ~] = fileparts(datFiles(i_dat_file).name);
        fileContents = fileread(fullfile(manually_identified_macro_contacts_path, datFiles(i_dat_file).name));
        lines = strsplit(fileContents, '\n');
        for row_id = 1:numel(lines)
            line = strtrim(lines{row_id});
            if contains(line, 'info') 
                break; % Skip this line
            end
            
            if ~isempty(line)
                coordinates = str2double(strsplit(line));
                x = coordinates(1);
                y = coordinates(2);
                z = coordinates(3);
                new_row = {shankID, row_id - skipped_line, x, y, z};
                contact_data = [contact_data; new_row];
            else
                skipped_line = skipped_line + 1;
            end
        end
    end
    contact_data_table = cell2table(contact_data, 'VariableNames', ...
        {'ShankID', 'ContactID', 'X_manual', 'Y_manual', 'Z_manual'});

elseif use_shank_model
    % get shank begin and end
    n_trajectory = length(imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory_id);
    all_traj_ids = string(imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory_id(:));
    n_shank = height(imaging_processing_info_struct.electrode_info);
    assert(n_trajectory == n_shank, ...
    'Mismatch: %d trajectories vs %d shanks in electrode_info.', n_trajectory, n_shank);
    contact_tbl_path = fullfile(imaging_processing_temp_folder, 'contact_data_table.mat');

    have_contact_tbl = false;
    if exist(contact_tbl_path, 'file')
        Sload = load(contact_tbl_path, 'contact_tbl');
        if isfield(Sload, 'contact_tbl') && istable(Sload.contact_tbl)
            contact_tbl = Sload.contact_tbl;
            have_contact_tbl = true;
        end
    end
    if ~have_contact_tbl
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
    end

    updated_imaging_processing_info_struct_save_path = [imaging_processing_temp_folder, ...
            '/imaging_processing_info_struct_after_electrode_localize.mat'];

    have_imaging_processing_info_struct_after_electrode_localize = false;
    if exist(updated_imaging_processing_info_struct_save_path, 'file')
        Sload = load(updated_imaging_processing_info_struct_save_path, 'imaging_processing_info_struct');
        if isfield(Sload, 'imaging_processing_info_struct') 
            tmp_imaging_processing_info_struct = Sload.imaging_processing_info_struct;
            have_imaging_processing_info_struct_after_electrode_localize = true;
            shank_model_all = tmp_imaging_processing_info_struct.shank_model_all;
        end
    end
    if ~have_imaging_processing_info_struct_after_electrode_localize
            shank_model_all = struct( ...
                'shank_id',{}, 'final_model',{}, 'delta_mm',{}, ...
                'proximal_ras',{}, 'distal_ras',{}, 'bezier_control',{}, ...
                'centerline_pts',{}, 'centerline_s',{}, ...
                'tube_radius_mm',{}, 'tube_auc_line',{}, 'tube_auc_bez',{}, ...
                'BIC_line',{}, 'BIC_bez',{}, 'contact_s_mm' ,{}, ...
                'contact_centers',{}, 'contact_axes',{}, 'contact_intensity',{},...
                'contact_score',{}, 'contact_center_CT_intensity',{}, 'contact_volume_mean_CT_intensity',{},...
                'contact_len_mm',{});
    end

    processed_from_models = strings(0,1);
    if ~isempty(shank_model_all)
        % shank_model_all is a struct array with field .shank_id
        processed_from_models = string({shank_model_all.shank_id}.');
    end
    processed_from_contacts = strings(0,1);
    if height(contact_tbl) > 0 && any(strcmpi(contact_tbl.Properties.VariableNames,'ShankID'))
        % Ensure column is string for set ops
        if ~isa(contact_tbl.ShankID,'string'); contact_tbl.ShankID = string(contact_tbl.ShankID); end
        processed_from_contacts = unique(contact_tbl.ShankID);
    end

    processed_from_models   = string({shank_model_all.shank_id});
    processed_from_contacts = string(contact_tbl.ShankID);
    
    % intersection: only shanks that are in BOTH sets
    processed_shanks = intersect(processed_from_models, processed_from_contacts);
    is_unprocessed = ~ismember(all_traj_ids, processed_shanks);
    unprocessed_traj_ids = all_traj_ids(is_unprocessed);
    fprintf('[Resume] %d/%d shanks already processed. %d remaining.\n', ...
        sum(~is_unprocessed), n_trajectory, sum(is_unprocessed));
    if isempty(unprocessed_traj_ids)
        warning('All shanks appear processed; nothing to do.');
        i_traj_list = [];   % nothing to loop
    else
        % Map unprocessed IDs back to indices in the original arrays
        i_traj_list = find(is_unprocessed);
    end

    for i_traj_unprocessed = 1:numel(i_traj_list)
        i_traj = i_traj_list(i_traj_unprocessed);
        traj_id = imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory_id{i_traj};
        % these are in mm, in RAS coordinates
        start_traj = imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory(i_traj).start;
        end_traj = imaging_processing_info_struct.coregistered_to_mri_trajectory.trajectory(i_traj).end;

        % viz
        figure('Color','w'); hold on;
        h = patch(isosurface(X, Y, Z, BW, 0.5));             % use meshgrid coords
        isonormals(X, Y, Z, double(ctVol), h);               % normals from intensity
        set(h, 'EdgeColor','none', 'FaceColor',[0.75 0.78 0.85], 'FaceAlpha',0.25);

        % OPTIONAL: decimate to speed up rendering
        % h = reducepatch(h, 0.5);  % keep 50% of faces
        % Transform vertices from voxel (1-based) -> 0-based -> RAS (mm)
        Vvox1 = h.Vertices;           % [i j k] in 1-based voxel indices
        Vras  = ([Vvox1, ones(size(Vvox1,1),1)] * Avox2ras1.');
        h.Vertices = Vras(:,1:3);

        % Scene aesthetics
        axis equal vis3d; camlight headlight; lighting gouraud;
        grid on; box on;
        xlabel('X_R_A_S (mm)'); ylabel('Y_R_A_S (mm)'); zlabel('Z_R_A_S (mm)');
        title(sprintf('CT ≥ %.0f (electrode isosurface in RAS)', thr));

        % plot the planned trajectory
        plot3([start_traj(1) end_traj(1)], [start_traj(2) end_traj(2)], [start_traj(3) end_traj(3)], ...
        '-', 'LineWidth', 2);

        % --- Endpoint spheres (approx macro size visually)
        [sx, sy, sz] = sphere(16);
        r = 0.6; % mm; tweak to ~0.6–1.0 to match your macro appearance
        surf(sx*r + start_traj(1), sy*r + start_traj(2), sz*r + start_traj(3), ...
            'EdgeColor','none','FaceAlpha',0.95);
        surf(sx*r + end_traj(1), sy*r + end_traj(2), sz*r + end_traj(3), ...
            'EdgeColor','none','FaceAlpha',0.95);

        % 3d rendering
        filename = [subj_id, '_ct_rotation_traj_', traj_id, '_and_metal_voxels'];
        full_path = fullfile(report_dir, [filename, '.mp4']);
        vObj = VideoWriter(full_path, 'MPEG-4');
        vObj.FrameRate = 30;       % frames per second
        open(vObj);

        nFrames = 360;             % one full circle (degrees)
        az0 = 135; el0 = 20;       % starting view angles

        for k = 1:nFrames
            az = az0 + (k-1);      % increment azimuth by 1° per frame
            view(az, el0);         % keep elevation constant, rotate around z
            drawnow;
            frame = getframe(gcf);
            writeVideo(vObj, frame);
        end

        close(vObj);
        fprintf('Saved rotation to %s\n', full_path);

        % get the entry point 
        % "proximal to distal" refers to the location
        % of the electrode contacts along the depth of the SEEG electrode
        % as it is inserted into the brain. "Proximal" refers to the
        % electrode contacts closer to the skull and brain surface, while
        % "distal" refers to those deeper within the brain tissue.

        % Sample along line in fine steps (e.g., 0.2 mm)
        % clear cache
        close all;

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
            warning('No skull intersection found on traj %s; check boneThr.', traj_id);
            entry_ras = start_traj; % use trajectory start point as fallback
        else
            entry_ras = P(hit,:); % proximal (on skull)
        end

        % diagnostic plot - entry pt here is at the planned trajectory
        % figure('Color','w'); hold on;
        % h = patch(isosurface(X, Y, Z, BW, 0.5));             % use meshgrid coords
        % isonormals(X, Y, Z, double(ctVol), h);               % normals from intensity
        % set(h, 'EdgeColor','none', 'FaceColor',[0.75 0.78 0.85], 'FaceAlpha',0.25);
        % 
        % % OPTIONAL: decimate to speed up rendering
        % % h = reducepatch(h, 0.5);  % keep 50% of faces
        % % Transform vertices from voxel (1-based) -> 0-based -> RAS (mm)
        % Vvox1 = h.Vertices;           % [i j k] in 1-based voxel indices
        % V0    = Vvox1 - 1;            % 0-based for vox2ras1
        % Vras  = ( [V0, ones(size(V0,1),1)] * Avox2ras.' );
        % h.Vertices = Vras(:,1:3);
        % 
        % % Scene aesthetics
        % axis equal vis3d; camlight headlight; lighting gouraud;
        % grid on; box on;
        % xlabel('X_R_A_S (mm)'); ylabel('Y_R_A_S (mm)'); zlabel('Z_R_A_S (mm)');
        % title(sprintf('CT ≥ %.0f (electrode isosurface in RAS)', thr));
        % 
        % % plot the planned trajectory
        % plot3([start_traj(1) end_traj(1)], [start_traj(2) end_traj(2)], [start_traj(3) end_traj(3)], ...
        % '-', 'LineWidth', 2);
        % 
        % % --- Endpoint spheres (approx macro size visually)
        % [sx, sy, sz] = sphere(16);
        % r = 0.6; % mm; tweak to ~0.6–1.0 to match your macro appearance
        % surf(sx*r + start_traj(1), sy*r + start_traj(2), sz*r + start_traj(3), ...
        %     'EdgeColor','none','FaceAlpha',0.95);
        % surf(sx*r + end_traj(1), sy*r + end_traj(2), sz*r + end_traj(3), ...
        %     'EdgeColor','none','FaceAlpha',0.95);
        % 
        % [sx, sy, sz] = sphere(40);
        % surf(sx*r + entry_ras(1), sy*r + entry_ras(2), sz*r + entry_ras(3), ...
        %     'EdgeColor','none','FaceAlpha',0.95);
        % 
        % plot3([start_traj(1) end_traj(1)], [start_traj(2) end_traj(2)], [start_traj(3) end_traj(3)], ...
        % '-', 'LineWidth', 2);


        roi_rad_mm = 5;
        ROI = cylinder_roi_mask(size(ctVol), Avox2ras0, start_traj, end_traj, roi_rad_mm);
        
        % Quick sanity plot: ROI + planned line (both in RAS)
        [m,n,p] = size(ctVol);
        [Xm,Ym,Zm] = meshgrid(1:n,1:m,1:p);
        figure('Color','w'); hold on;
        h_roi = patch(isosurface(Xm, Ym, Zm, ROI, 0.5)); set(h_roi,'EdgeColor','none','FaceAlpha',0.25);
        % -> RAS
        V0 = h_roi.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; h_roi.Vertices = Vras(:,1:3);
        plot3([start_traj(1) end_traj(1)], [start_traj(2) end_traj(2)], [start_traj(3) end_traj(3)], 'k--','LineWidth',2);
        axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
        xlabel('X_RAS (mm)'); ylabel('Y_RAS (mm)'); zlabel('Z_RAS (mm)');
        title('Exact cylinder ROI + planned trajectory'); view(135,20);
        filename = [subj_id, '_cylinder_around_planned_traj_lead_', traj_id];
        full_path = fullfile(report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');

        % sanity check on metal + ROI
        metal_roi = BW & ROI;
        figure('Color','w'); hold on;
        h_mroi = patch(isosurface(Xm, Ym, Zm, metal_roi, 0.5));
        set(h_mroi,'EdgeColor','none','FaceColor',[0.9 0.2 0.2],'FaceAlpha',0.45);
        V0 = h_mroi.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; h_mroi.Vertices = Vras(:,1:3);
        plot3([start_traj(1) end_traj(1)], [start_traj(2) end_traj(2)], [start_traj(3) end_traj(3)], 'k--','LineWidth',2);
        axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
        xlabel('X_RAS'); ylabel('Y_RAS'); zlabel('Z_RAS'); view(135,20);
        title('Metal ∩ ROI + planned trajectory');
        xlim([-80, 80]);
        ylim([-90, 90]);
        zlim([-90, 90]);
        filename = [subj_id, '_isolate_metal_close_to_plan_trajectory_lead_', traj_id];
        full_path = fullfile(report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');

        % combine plot
        figure('Color','w'); hold on;
        h_metal = patch(isosurface(Xm, Ym, Zm, BW, 0.5));
        isonormals(Xm, Ym, Zm, double(ctVol), h_metal);
        set(h_metal, 'EdgeColor','none', 'FaceColor', [0.80 0.20 0.20], 'FaceAlpha', 0.20);
        V0 = h_metal.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; h_metal.Vertices = Vras(:,1:3);

        % Metal and ROI surface (highlighted)
        metal_roi = BW & ROI;
        h_roi = patch(isosurface(Xm, Ym, Zm, metal_roi, 0.5));
        set(h_roi, 'EdgeColor','none', 'FaceColor', [0.10 0.45 0.85], 'FaceAlpha', 0.45);
        V0 = h_roi.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; h_roi.Vertices = Vras(:,1:3);
        plot3([start_traj(1) end_traj(1)], [start_traj(2) end_traj(2)], [start_traj(3) end_traj(3)], 'k--','LineWidth',2);
        axis equal vis3d; camlight headlight; lighting gouraud;
        grid on; box on;
        xlabel('X_RAS (mm)'); ylabel('Y_RAS (mm)'); zlabel('Z_RAS (mm)');
        title('Metal mask (red), Metal∩ROI (blue), Planned trajectory (black)');
        legend([h_metal h_roi], {'Metal BW','Metal ∩ ROI'}, 'Location','northeast');
        view(135,20);
        filename = [subj_id, '_isolate_metal_close_to_plan_trajectory_as_part_of_metal_voxels_lead_', traj_id];
        full_path = fullfile(report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');


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
        topK         = 5;          % how many seeds to keep (try 4–10)
        min_gap_mm   = 4.0;        % enforce spacing between seeds in mm
        sigma_vox    = 0.6;        % smoothing for peak picking (in voxels)
        dmax_mm = roi_rad_mm / 1.5;
        wI          = 0.7;         % weight: brightness
        wD          = 0.3;         % weight: distance (closer is better)

        
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
        
        % Normalize intensity & distance to [0,1]
        if isempty(vals)
            seed_ras = zeros(0,3); seed_ijk = zeros(0,3); seed_vals = [];
        else
            vmin = prctile(vals, 5); vmax = prctile(vals, 95);
            Inorm = (vals - vmin) / max(vmax - vmin, eps); Inorm = min(max(Inorm,0),1);
            Dnorm = 1 - min(rad_mm / dmax_mm, 1);   % 1 near line, 0 at dmax
        
            % Composite score: higher is better
            score = wI*Inorm + wD*Dnorm;
        
            % Sort by score desc
            [~, ord] = sort(score, 'descend');
            iPk = iPk(ord); jPk = jPk(ord); kPk = kPk(ord);
            cand_ras = cand_ras(ord,:); score = score(ord); vals = vals(ord);
        
            % Greedy non‑max suppression in **mm** space (min_gap_mm)
            chosen = false(size(iPk));
            seed_ras = zeros(0,3); seed_vals = [];
            for t = 1:numel(iPk)
                Pr = cand_ras(t,:);
                if isempty(seed_ras)
                    accept = true;
                else
                    d = sqrt(sum((seed_ras - Pr).^2,2));
                    accept = all(d >= min_gap_mm);
                end
                if accept
                    chosen(t)  = true;
                    seed_ras   = [seed_ras; Pr];         %#ok<AGROW>
                    seed_vals  = [seed_vals; vals(t)];    %#ok<AGROW>
                    if size(seed_ras,1) >= topK, break; end
                end
            end
        
            iSeed = iPk(chosen); jSeed = jPk(chosen); kSeed = kPk(chosen);
            seed_ijk = [iSeed, jSeed, kSeed];             % 1-based voxel indices
        end

        figure('Color','w'); hold on;

        [m,n,p] = size(ctVol);
        [Xm, Ym, Zm] = meshgrid(1:n, 1:m, 1:p);
        S = isosurface(Xm, Ym, Zm, metal_roi_lo, 0.5);      % voxel(1-based) space verts
        
        % Map ROI vertices to RAS using 0-based affine
        V0   = S.vertices - 1;                               % 1-based -> 0-based
        Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.';       % -> RAS (mm)
        S.vertices = Vras(:,1:3);
        
        % --- Sample CT intensity at those vertices for coloring ---
        % Convert RAS verts -> voxel(1-based, continuous) via vox2ras1
        ijk1_h = (Avox2ras1 \ [S.vertices, ones(size(S.vertices,1),1)].').'; % Nx4
        ijk1   = ijk1_h(:,1:3);  % [i j k] (row,col,slice), continuous
        
        % Trilinear sample
        I = trilinear(double(ctVol), ijk1);                  % size = [Nverts x 1]
        
        % Normalize to [0,1] for color & alpha
        Imin = prctile(I, 5); Imax = prctile(I, 99);
        C = (I - Imin) / max(Imax - Imin, eps);
        C = min(max(C,0),1);
        
        % Make alpha a softer version of color (tweak factor if you like)
        A = ones(size(C)) * 0.2;                                  % between ~0.15 and 0.70
        

        h_roi = patch('Vertices', S.vertices, 'Faces', S.faces, ...
                      'FaceVertexCData', C, ...
                      'FaceColor','interp', ...              % color by intensity
                      'EdgeColor','none', ...
                      'FaceVertexAlphaData', A, ...
                      'FaceAlpha','interp');                 % transparency by intensity
        colormap turbo; colorbar; clim([0 1]);
        
        % Planned line (reference)
        plot3([start_traj(1) end_traj(1)], ...
              [start_traj(2) end_traj(2)], ...
              [start_traj(3) end_traj(3)], 'k--','LineWidth',1.2);
        
        % --- Seeds: colored by composite score ---
        % Assumes you already computed:
        %   seed_ras  (Ns x 3, RAS mm)
        %   score     (Nc x 1 for all candidates)
        %   chosen    (Nc x 1 logical indicating which candidates became seeds)
        % If you only have seed_ras and per-seed composite scores, just use them.
        
        % Plot all candidate peaks faintly (optional)
        % scatter3(cand_ras(:,1), cand_ras(:,2), cand_ras(:,3), 10, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha',0.25);
        
        % Plot chosen seeds colored by their composite score
        scatter3(seed_ras(:,1), seed_ras(:,2), seed_ras(:,3), 64, score(chosen), 'filled', ...
                 'MarkerEdgeColor','k', 'LineWidth',0.8);
        % Keep same turbo colormap; rescale colorbar to seed score range if you want:
        % caxis([min(score(chosen)) max(score(chosen))]);
        
        axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
        xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
        title('Metal ROI colored by CT intensity; seeds colored by composite score');
        view(135,20);
        filename = [subj_id, '_get_seeds_for_growing_to_identify_lead_', traj_id];
        full_path = fullfile(report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');

        % growing algorithm
        % we have multiple seeds, we will connect the neighbors,
        % to only connect voxels corresponding a lead, we will limit
        % growing direction. In practice, we will continue to grow, and
        % along the process, we will get a tube following the planned
        % trajectory direction that encloses all the connected voxels.
        % we can somehow assign a score for each voxels to represent the
        % percentage of filling the tube, if belong to the same lead, the
        % filling percentage is high while if belonging to different lead,
        % the percentage will be low. we should start the calculation of
        % filling percentage froom 2 voxels.In this way, we probably do not
        % need bridge, instead of we will search for closest point from the
        % growing bolb to decide what is next to include


        [M,N,P] = size(ctVol);
        [Xras, t_field] = precompute_ras_t(Avox2ras0, start_traj, U, M,N,P);
        
        vote_frac   = 0.30;     % fraction of seeds that must vote for a voxel, at least included in 2 volumes
        min_cc_vox  = 50;       % drop tiny components after consensus
        max_voxels  = 2e6;
        
        % grow_from_seed options (adjust as you like)
        opt.UseROI      = false;
        opt.AngleDeg    = 3;
        opt.TPadMM      = 5.0;
        opt.Weights     = [1 0.1];   % [alpha_fill, beta_dist]
        opt.ScoreMin    = -0.3;
        opt.fill_drop_tol = 0.2;
        opt.FillMin     = 0.05;
        opt.KNearest    = 1;
        opt.KInit       = 8;
        opt.InitUseROI  = false;
        
        % Voxel sizes [I J K] in mm
        vx = [imaging_processing_info_struct.ct_volume_registered_to_mri.ysize, ...
              imaging_processing_info_struct.ct_volume_registered_to_mri.xsize, ...
              imaging_processing_info_struct.ct_volume_registered_to_mri.zsize];
        
        % --------- Parallel growth over seeds ----------
        K = size(seed_ijk, 1);
        Glist     = cell(K,1);
        info_list = cell(K,1);
        
        % 
        parfor ss = 1:K
            seed = seed_ijk(ss,:);

            [Gss, infoss] = grow_from_seed( ...
                seed, ctVol, ROI, metal_roi_lo, Avox2ras0, ...
                start_traj, U, max_voxels, vx, ...
                'UseROI',      opt.UseROI, ...
                'AngleDeg',    opt.AngleDeg, ...
                'TPadMM',      opt.TPadMM, ...
                'Weights',     opt.Weights, ...
                'ScoreMin',    opt.ScoreMin, ...
                'FillMin',     opt.FillMin, ...
                'KNearest',    opt.KNearest, ...
                'KInit',       opt.KInit, ...
                'InitUseROI',  opt.InitUseROI);
        
            Glist{ss}     = logical(Gss);
            info_list{ss} = infoss;
        end
        
        % consensus voting
        vote = zeros(size(ctVol), 'uint16');
        for ss = 1:K
            if ~isempty(Glist{ss})
                vote = vote + uint16(Glist{ss});
            end
        end
        
        maj_thr = ceil(vote_frac * K);     
        metal_roi_unique = vote >= maj_thr;

        % add a assertion here, if more than 
        assert( (nnz(metal_roi_unique) / nnz(metal_roi)) > 0.5, ['identified lead is less than half' ...
            ' of the size of metal surrounding planned traj'])

                

        % remove tiny islands
        % metal_vote = bwareaopen(metal_vote, min_cc_vox);
        % 
        % % ---------- pick dominant connected component (avoid merging another lead) ----------
        % CC = bwconncomp(metal_vote, 2);
        % if CC.NumObjects >= 2
        %     % score each CC by (overlap with high-threshold metal) and axial span along U
        %     scores = zeros(CC.NumObjects,1);
        %     spans  = zeros(CC.NumObjects,1);
        %     for c = 1:CC.NumObjects
        %         idx = CC.PixelIdxList{c};
        %         overlap = sum(metal_roi_hi(idx));  % favors truly metal-dense component
        %         [ii,jj,kk] = ind2sub(size(metal_vote), idx);
        %         % project to line (t in mm) using voxel centers (0-based -> RAS)
        %         P = (Avox2ras0 * [jj-1, ii-1, kk-1, ones(numel(idx),1)].').';
        %         P = P(:,1:3);
        %         t = dot(P - start_traj, repmat(U,size(P,1),1), 2);
        %         spans(c)  = prctile(t,95) - prctile(t,5);
        %         scores(c) = 1.0*overlap + 0.02*spans(c);  % weights: tweak if needed
        %     end
        %     [~,bestC] = max(scores);
        %     metal_roi_unique = false(size(metal_vote));
        %     metal_roi_unique(CC.PixelIdxList{bestC}) = true;
        % else
        %     metal_roi_unique = metal_vote;
        % end

        
        % (optional) light cleanup confined to tube
        % tube_mask = cylinder_roi_mask(size(ctVol), Avox2ras0, start_traj, end_traj, r_gate_mm + 0.5);
        % metal_roi_unique = metal_roi_unique & tube_mask;
        % metal_roi_unique = bwareaopen(metal_roi_unique, 50);  % drop tiny crumbs


        % viz
        % combine plot
        figure('Color','w'); hold on;
        h_metal = patch(isosurface(Xm, Ym, Zm, BW, 0.5));
        % make 3d surface look smooth
        isonormals(Xm, Ym, Zm, double(ctVol), h_metal);
        set(h_metal, 'EdgeColor','none', 'FaceColor', [0.80 0.20 0.20], 'FaceAlpha', 0.20);
        V0 = h_metal.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; h_metal.Vertices = Vras(:,1:3);

        h_roi = patch(isosurface(Xm, Ym, Zm, metal_roi_unique, 0.5));
        set(h_roi, 'EdgeColor','none', 'FaceColor', [0.10 0.45 0.85], 'FaceAlpha', 0.45);
        V0 = h_roi.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; h_roi.Vertices = Vras(:,1:3);
        plot3([start_traj(1) end_traj(1)], [start_traj(2) end_traj(2)], [start_traj(3) end_traj(3)], 'k--','LineWidth',2);
        axis equal vis3d; camlight headlight; lighting gouraud;
        grid on; box on;
        xlabel('X_RAS (mm)'); ylabel('Y_RAS (mm)'); zlabel('Z_RAS (mm)');
        title('Metal mask (red), Metal∩ROI (blue), Planned trajectory (black)');
        legend([h_metal h_roi], {'Metal BW','Lead ∩ ROI'}, 'Location','northeast');
        view(135,20);
        filename = [subj_id, '_isolate_lead_', traj_id, '_from_metal_close_to_plan_trajectory_as_part_of_metal_voxels'];
        full_path = fullfile(report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');

        % now we identify metal volume corresponding to the shank
        % we will define shank entry and end pt, and potentialy the third
        % point, the structure will be added to     
        % shank_models = struct('shank_id',{}, 'proximal_ras',{}, 
        % 'distal_ras',{}, 'dir_unit',{}, 'third_flex_pt',{});
        % we will start with linear fitting as a prior, than do a 3
        % point-intepolation spilne see if it fits the vlume better
        %, approach for quantifying better could be 'baysien factor',
        % the start pt is the one close to the skull
        bm_path = [data_path, subj_id, '/IMAGING/segmentation/mri/brainmask.auto.mgz'];
        % Load brainmask
        bm = MRIread(bm_path);
        out = fit_shank_line_from_blob(metal_roi_unique, BW_skull, Avox2ras0, start_traj, end_traj, ...
            'BrainMaskMGZ', bm, 'BrainMaskVox2RAS', bm.vox2ras0, 'VisualizeBands', true, ...
            'subj_id', subj_id, 'report_dir', report_dir, 'traj_id', traj_id);
        % when fitting the bolb, we only consider the part within the skull
        % below code is for diagnostic 
        
        % bm_path = [data_path, subj_id, '/IMAGING/segmentation/mri/brainmask.auto.mgz'];
        % % Load brainmask
        % bm = MRIread(bm_path);
        % BW_brain = bm.vol > 0;
        % 
        % A_mri = bm.vox2ras;      % MRI vox->RAS
        % 
        % sz_ct = size(metal_roi_unique);       % CT grid size
        % 
        % BW_brain_in_CT = resample_mask_nn(BW_brain, A_mri, sz_ct, imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras0);
        % 
        % metal_roi_unique_within_skull = metal_roi_unique & BW_brain_in_CT;
        % metal_roi_unique_outside_skull = metal_roi_unique & (~metal_roi_unique_within_skull);
        % 
        % fprintf('Kept %d voxels inside brain mask (%.1f%% of metal)\n', ...
        %     nnz(metal_roi_unique_within_skull), ...
        %     100*nnz(metal_roi_unique_within_skull)/nnz(metal_roi_unique));
        % 
        % 
        % [m,n,p] = size(metal_roi_unique);
        % [Xm, Ym, Zm] = meshgrid(1:n, 1:m, 1:p);
        % toRAS = @(V1) ([V1-1, ones(size(V1,1),1)] * Avox2ras.');
        % 
        % figure('Color','w'); hold on;
        % 
        % % Skull
        % if any(BW(:))
        %     hs = patch(isosurface(Xm, Ym, Zm, BW, 0.5));
        %     set(hs,'EdgeColor','none','FaceColor',[0.85 0.80 0.70],'FaceAlpha',0.30);
        %     V0 = hs.Vertices; Vras = toRAS(V0); hs.Vertices = Vras(:,1:3);
        % end
        % 
        % % Periosteal band metal (outside skull layer)
        % if any(metal_roi_unique_within_skull(:))
        %     hperi = patch(isosurface(Xm, Ym, Zm, metal_roi_unique_within_skull, 0.5));
        %     set(hperi,'EdgeColor','none','FaceAlpha',0.65,'FaceColor',[0.90 0.20 0.20]); % red
        %     V0 = hperi.Vertices; Vras = toRAS(V0); hperi.Vertices = Vras(:,1:3);
        % end
        % 
        % % Cortical-outer band metal (within skull outer layer)
        % if any(metal_roi_unique_outside_skull(:))
        %     hcort = patch(isosurface(Xm, Ym, Zm, metal_roi_unique_outside_skull, 0.5));
        %     set(hcort,'EdgeColor','none','FaceAlpha',0.60,'FaceColor',[0.10 0.45 0.85]); % blue
        %     V0 = hcort.Vertices; Vras = toRAS(V0); hcort.Vertices = Vras(:,1:3);
        % end
        % 
        % 
        % axis equal vis3d; camlight headlight; lighting gouraud; grid on; box on;
        % xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
        % title(sprintf('Skull (tan), Periosteal metal (red, ≤%.1f mm), Cortical-outer metal (blue, ≤%.1f mm)', band_mm, band_mm));

        % electrode identification
        % the idea would be find a predefined electrode (we know their relative location)
        % for DIXI, each elecrtode length is associated with a configuration
        % for example: elec_model_dixi_15cm = [0 3.5 7 10.5 14 27 30.5 34 37.5 41 54 57.5 61 64.5 68];
        % elec_model_dixi_18cm = [0 3.5 7 10.5 14 17.5 30.5 34 37.5 41 44.5 48 61 64.5 68 71.5 75 78.5];
        % elec_model_dixi_am = 0:3.5:59.5;
       

        % Isotropic resampling - opt out because it introduced disaligned
        % voxels
        % iso_mm = min(round(vx, 1));
        % 
        % dx = norm(Avox2ras0(1:3,1));
        % dy = norm(Avox2ras0(1:3,2));
        % dz = norm(Avox2ras0(1:3,3));
        % 
        % sz = size(ctVol);
        % 
        % % --- 8 corners of the source volume in voxel coords (0-based) ---
        % corn = [
        %     0       0       0
        %     sz(1)-1 0       0
        %     0       sz(2)-1 0
        %     0       0       sz(3)-1
        %     sz(1)-1 sz(2)-1 0
        %     sz(1)-1 0       sz(3)-1
        %     0       sz(2)-1 sz(3)-1
        %     sz(1)-1 sz(2)-1 sz(3)-1
        % ];
        % corn_h = [corn ones(8,1)];
        % RAScorn = (Avox2ras0.' \ corn_h.').';  % [x y z 1] = [u v w 1] * inv(A)^T
        % RAScorn = RAScorn(:,1:3);
        % 
        % % --- Tight RAS-aligned box around the oblique volume ---
        % xmin = -150; xmax = 150;
        % ymin = -150; ymax = 100;
        % zmin = -150; zmax = 150;
        % 
        % xg = xmin:iso_mm:xmax;
        % yg = ymin:iso_mm:ymax;
        % zg = zmin:iso_mm:zmax;
        % [Xras, Yras, Zras] = ndgrid(xg, yg, zg);  % sizes: [numel(xg) x numel(yg) x numel(zg)]
        % 
        % % 2) Map RAS -> source voxel (0-based)
        % AinvT = inv(Avox2ras0).';                       % for row-vector multiply
        % UVW   = [Xras(:) Yras(:) Zras(:) ones(numel(Xras),1)] * AinvT;  % -> [i j k 1]
        % I0 = reshape(UVW(:,1), size(Xras));   % 0-based i
        % J0 = reshape(UVW(:,2), size(Xras));   % 0-based j
        % K0 = reshape(UVW(:,3), size(Xras));   % 0-based k
        % 
        % % 3) Interpolate (remember: F expects (J,I,K) in 1-based coords)
        % F = griddedInterpolant(double(ctVol), 'linear', 'none');
        % ctVolIso = F(J0+1, I0+1, K0+1);       % <-- this ordering is crucial
        % % ctVolIso = F(I0+1, J0+1, K0+1); 
        % 
        % % 4) New affine for the regular RAS grid (diagonal spacing, origin at first samples)
        % Avox2ras0Iso = [iso_mm 0      0      xg(1);
        %                 0      iso_mm 0      yg(1);
        %                 0      0      iso_mm zg(1);
        %                 0      0      0      1];

        % validation
        % [i1,j1,k1] = deal(0,0,0);
        % pIso = Avox2ras0Iso * [i1;j1;k1;1];  % should be [xg(1); yg(1); zg(1); 1]
        % assert(all(abs(pIso(1:3) - [xg(1); yg(1); zg(1)]) < 1e-9));
        % 
        % [i2,j2,k2] = deal(5,7,9);
        % pIso = Avox2ras0Iso * [i2;j2;k2;1];
        % assert(all(abs(pIso(1:3) - [xg(i2+1); yg(j2+1); zg(k2+1)]) < 1e-9), 'Avox2ras0Iso mismatch');

        % output to check - verified in freesurfer
        % fn_iso = fullfile(imaging_processing_temp_folder, sprintf('%s_ct_iso_%0.1fmm.nii.gz', subj_id, iso_mm));
        % 
        % % 1) Build a NIfTI object
        % dat = single(ctVolIso);                       % or uint16(max(0,min(65535,round(ctVolIso))))
        % nii = make_nii(dat, [iso_mm iso_mm iso_mm], [0 0 0], 16);  % 16=float32; use 512 for uint16
        % 
        % nii.hdr.hist.sform_code = 1;   % scanner/anatomical
        % nii.hdr.hist.qform_code = 0;   % disable qform so sform is authoritative
        % 
        % % 3) Copy your affine rows into srow_x/y/z
        % % Avox2ras0Iso is:
        % %   [ iso  0    0    x0
        % %     0    iso  0    y0
        % %     0    0    iso  z0
        % %     0    0    0    1 ]
        % nii.hdr.hist.srow_x = Avox2ras0Iso(1,:);
        % nii.hdr.hist.srow_y = Avox2ras0Iso(2,:);
        % nii.hdr.hist.srow_z = Avox2ras0Iso(3,:);
        % 
        % nii.hdr.hist.intent_name = sprintf('Isotropic CT %.1f mm (RAS)', iso_mm);
        % 
        % % 4) Save
        % save_nii(nii, fn_iso);
        

        if electrode_manufacture == 'DIXI'
            contact_radius_mm = DIXI_am_diameter/2;
            n_contact_str = imaging_processing_info_struct.electrode_info.n_contact_str{i_traj}; % if exists
            n_contact_num = imaging_processing_info_struct.electrode_info.n_contact(i_traj);     % if exists
            mdl = dixi_model_from_row(n_contact_str, n_contact_num);
        end

        if strcmpi(out.model, 'bezier') && ~isempty(out.bezier_control)
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
        pt_at_s = @(ss) [interp1(s, Cpts(:,1), ss, 'linear','extrap'), ...
                         interp1(s, Cpts(:,2), ss, 'linear','extrap'), ...
                         interp1(s, Cpts(:,3), ss, 'linear','extrap')];

        rad_mm = mdl.diameter_mm/2; 
        len_mm = mdl.contact_len_mm;
        offs   = mdl.offsets_mm(:)';        % from distal tip toward proximal (mm)
        
        % “s” location of distal tip on the centerline:
        % By construction above, s=0 at proximal (E) and s=Ltot at distal (D).
        s_distal = Ltot;
         % distal tip (deepest point inside brain).

        % robust method
        score_pose = @(C,u) score_pose_contrast(C,u,rad_mm,len_mm,Avox2ras0,double(ctVol),vx);
        ds = min(vx)/2;                                % mm step
        s_grid = (0:ds:Ltot + DIXI_am_spacing/2).';
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

        % 1) Global shift δ (mm) to align distal offsets to the score profile

        % Intersect with your desired search range
        delta_lo = -DIXI_am_spacing/2;
        delta_hi = DIXI_am_spacing/2;
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
        p_of_s_safe = @(ss) deal_vals(ss, F, 0, Ltot + DIXI_am_spacing);                  % -Inf outside [0,Ltot]
        
        % Expected contact centers given a shift dlt
        exp_centers = @(dlt) s_distal - (offs(:) + len_mm/2 - dlt);     % column vector
        
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
        
        figure('Color','w'); tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
        
        % (1) Objective vs shift
        nexttile;
        plot(delta_grid, obj_grid_orig, '-', 'LineWidth',1.2); hold on;
        plot(delta_grid, obj_grid_rob,  '--', 'LineWidth',1.4);
        xline(best_delta_orig, ':', sprintf('\\delta_{orig}=%.3f', best_delta_orig), 'LabelOrientation','aligned');
        xline(best_delta_rob , '-.', sprintf('\\delta_{rob}=%.3f',  best_delta_rob ), 'LabelOrientation','aligned');
        xlabel('\delta (mm)'); ylabel('Objective');
        legend('Original sum','Trimmed (robust)','Location','best'); grid on;
        title('Global shift fit (objective vs. \delta)');
        
        % (2) p_det(s) with expected centers from both fits
        nexttile;
        plot(s_grid, p_det, 'k-', 'LineWidth',1.0); hold on;
        yl = ylim;
        stem_safe(s_exp_orig, p_of_s_safe(s_exp_orig), 'Marker','o','Color',[0.10 0.45 0.85], 'DisplayName','centers (orig)');
        stem_safe(s_exp_rob , p_of_s_safe(s_exp_rob ), 'Marker','^','Color',[0.85 0.25 0.10], 'DisplayName','centers (rob)');
        ylim(yl);
        xlabel('s (mm)'); ylabel('p_{det}(s)');
        legend('Location','best'); grid on;
        title('Detrended profile with expected centers');

        delta_use = best_delta_orig;      % or: best_delta_rob

        % Expected centers from the chosen global delta (comb)
        s_seed = exp_centers(delta_use);                 % N×1
        s_seed = min(max(s_seed,0), Ltot + DIXI_am_spacing / 2);               % clip to domain
        
        % --- PER-CONTACT MICRO REFINE (around each seed)
        ds     = s_grid(2)-s_grid(1);
        ax_win = DIXI_am_spacing / 20; 
        
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
        nexttile(2); % second panel in the tiledlayout
        s_refined_orig = s_exp_orig;
        stem_safe(s_refined, p_of_s_safe(s_refined), 'Marker','o','Color',[0.05 0.90 0], 'DisplayName','refined (orig)');



        contact_s = s_refined_orig(:);                           % Nx1

        cellsC = arrayfun(@(ss) pt_at_s(ss), contact_s, 'UniformOutput', false);
        contact_centers = cat(1, cellsC{:});     % or: vertcat(cellsC{:})
        
        % Tangents (N×3)
        cellsU = arrayfun(@(ss) tan_at_s(ss, s, Cpts), contact_s, 'UniformOutput', false);
        contact_tangents = cat(1, cellsU{:});
        
        % Along-track score (N×1)
        contact_score = p_of_s_safe(contact_s);
        
        % CT intensities (nearest voxel + small orthogonal disc average),
        % in sample_ct_nn we must do ijk index
        contact_HU_nn   = sample_ct_nn(double(ctVol), Avox2ras0, contact_centers);           % N×1
        contact_HU_disc = arrayfun(@(k) meanHU_disc(double(ctVol), Avox2ras0, ...
                                       contact_centers(k,:), contact_tangents(k,:), ...
                                       0.6*(mdl.diameter_mm/2), 0.6*len_mm), (1:numel(contact_s))');



        filename = [subj_id, '_identify_contact_from_lead_', traj_id, '_anulus_score'];
        full_path = fullfile(report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');

        % note Vq = interp3(V,Xq,Yq,Zq) assumes a default grid of sample
        % points. The default grid points cover the region, X=1:n, Y=1:m,
        % Z=1:p, where [m,n,p] = size(V). Use this syntax when you want to
        % conserve memory and are not concerned about the absolute
        % distances between points. we have verified cpts are correct RAS
        % coordinates meshgrid's inversion of the first two indices can
        % help when dealing with plotting, or with some of the output of
        % Image Processing Toolbox functions.

        % However, when you want to calculate, say, the values of a 4D Gaussian,
        % your only choice is to use ndgrid:

        % if doing ndgrid, it returns i(row) j(col) k, later when we do isosurface,
        % it needs  col row ..
        if exist('metal_roi_unique','var') && any(metal_roi_unique(:))
            [m,n,pz] = size(metal_roi_unique);
            [I,J,K] = meshgrid(1:m, 1:n, 1:pz);
            figure('Color','w'); hold on; axis equal vis3d; xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Metal isosurface with contact centers');
            hs = patch(isosurface(I, J, K, metal_roi_unique, 0.2));
            set(hs,'EdgeColor','none','FaceAlpha',0.5);
            Vijk = hs.Vertices;
            Vras = ([Vijk - 1, ones(size(Vijk,1),1)] * Avox2ras0.');   % -> [x y z 1]
            hs.Vertices = Vras(:,1:3);
            hline = plot3(Cpts(:,1), Cpts(:,2), Cpts(:,3), ...
                'w-', 'LineWidth', 0.8, 'DisplayName','centerline', 'Color', ...
                'green');
            uistack(hline,'top');  % ensure line is on top
            % --- Direct voxel HU at nearest voxels to the isosurface vertices ---
            % Vijk(:,1) is j (column)
            Ii = round(Vijk(:,2));  Jj = round(Vijk(:,1));  Kk = round(Vijk(:,3));
            Ii = max(1,min(m,Ii));  Jj = max(1,min(n,Jj));  Kk = max(1,min(pz,Kk));
            lin_idx   = sub2ind([m n pz], Ii, Jj, Kk);
            HU_direct = double(ctVol(lin_idx));
        
            % Colorize by CT HU
            % interp3 uses subscripts as I J K?
            Ainv = inv(Avox2ras0);
            uvw1 = [Vras(:,1:3), ones(size(Vras,1),1)] * Ainv.';       % -> 0-based voxel
            Iq = uvw1(:,1) + 1;  Jq = uvw1(:,2) + 1;  Kq = uvw1(:,3) + 1;
            CtV = interp3(double(ctVol), Iq, Jq, Kq, 'linear', NaN);
            % figure; scatter(HU_direct, CtV);
            % If we do scatter, they are correlated.
            % If we use this CtV, we will only show the intensity of the
            % surface, later it is more like the intensity of the actual
            % voxels

            bad = isnan(CtV);
            if any(bad)
                CtV(bad) = interp3(double(ctVol), Iq(bad), Jq(bad), Kq(bad), 'nearest', 0);
            end
            % Precompute nearest-metal index map once per figure
            % [~, idxNearestMetal] = bwdist(metal_roi_unique);    % idx of nearest nonzero voxel
            % 
            % % Round the vertex sample point to the nearest voxel
            % Ii = round(Iq);  Jj = round(Jq);  Kk = round(Kq);
            % Ii = max(1, min(m, Ii));  Jj = max(1, min(n, Jj));  Kk = max(1, min(pz, Kk));
            % lin = sub2ind([m n pz], Ii, Jj, Kk);
            % 
            % % Jump to the nearest metal voxel and sample HU there
            % linMetal = idxNearestMetal(lin);                    % nearest voxel inside metal blob
            % CtV = double(ctVol(linMetal));

            hs.FaceVertexCData = CtV;
            hs.FaceColor = 'interp'; colormap(parula); colorbar;
            hu_lo = prctile(CtV, 60);    % below this → very transparent
            hu_hi = prctile(CtV, 90);   % at/above this → almost opaque
            alpha = (CtV - hu_lo) ./ max(1, (hu_hi - hu_lo));   % normalize
            alpha = min(max(alpha, 0), 1);                      % clamp [0,1]
            alpha = 0.20 + 0.70*alpha;                          % map to [0.10, 0.90]
            
            hs.FaceVertexAlphaData = alpha;
            hs.FaceAlpha = 'interp';          % interpolate vertex alphas across faces
            hs.AlphaDataMapping = 'none';
            clim([hu_lo hu_hi]); 
                    
            % Compute RAS coordinates of centers
            cells  = arrayfun(@(ss) pt_at_s(ss), s_refined_orig(:), 'uni', 0);
            C_orig = vertcat(cells{:}); 
            % C_rob  = cell2mat(arrayfun(@(ss) pt_at_s(ss), s_refined_rob , 'uni', 0).');
        
            % Plot spheres (or markers) for both sets
            plot_contact_markers_3d(C_orig, 0.3, [1, 0, 0], 'orig');  % blue-ish
            % plot_contact_markers_3d(C_rob , 0.6, [0.85 0.25 0.10], 'rob');   % red-ish
        
            camlight headlight; lighting gouraud; grid on; box on;

            % --- draw proximal (entry) and distal points ---
            % plot3(out.proximal_entry_ras(1), out.proximal_entry_ras(2), out.proximal_entry_ras(3), ...
            %       'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName','proximal entry');
            % 
            % plot3(out.distal_ras(1), out.distal_ras(2), out.distal_ras(3), ...
            %       'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName','distal');

            
            % plot3(Cpts(end,1), Cpts(end,2), Cpts(end,3), ...
            %       'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName','distal');

            % save figure here
            filename = [subj_id, '_identify_contact_from_lead_', traj_id, '_colorcode_surface_intensity'];
            full_path = fullfile(report_dir, [filename, '.pdf']);
            saveas(gcf, full_path, 'pdf');
        end

        % verify if the transformation is correct
        % [ix, iy, iz] = ind2sub(size(metal_roi_unique), find(metal_roi_unique, 1));
        % 
        % % Original HU directly from ctVol
        % HU_direct = double(ctVol(ix, iy, iz));
        % 
        % % Convert voxel(1-based) -> voxel(0-based) coords
        % vox0 = [ix-1, iy-1, iz-1, 1];             % row vector
        % ras  = vox0 * Avox2ras0.';                % [x y z 1] in RAS
        % 
        % % Now inverse: RAS -> voxel(0-based)
        % uvw0 = [ras(1:3), 1] * inv(Avox2ras0).';  % back to [i j k 1] in 0-based
        % ii = uvw0(1)+1; jj = uvw0(2)+1; kk = uvw0(3)+1;  % convert to 1-based for interp3
        % 
        % % Interpolated HU at returned voxel coords
        % HU_interp = interp3(double(ctVol), jj, ii, kk, 'linear', NaN);
        % fprintf('Direct HU = %.1f, round-trip interp HU = %.1f\n', HU_direct, HU_interp);


        % clear cache
        close all;

        % diagnostic, - metal blob with
        figure('Color','w'); hold on;

        h_roi = patch(isosurface(Xm, Ym, Zm, metal_roi_unique, 0.5));
        set(h_roi, 'EdgeColor','none', 'FaceColor', [0.10 0.45 0.85], 'FaceAlpha', 0.2);
        V0 = h_roi.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; h_roi.Vertices = Vras(:,1:3);
        nTheta = 28; nSeg = 6;   % rendering resolution (tweak for speed)
        for k = 1:size(contact_centers,1)
            Ck = contact_centers(k,:); uk = contact_tangents(k,:);
            [V,F] = cylinder_at_pose(Ck, uk, len_mm, rad_mm, nTheta, nSeg);
            if ~isempty(V)
                patch('Vertices',V,'Faces',F, ...
                    'FaceColor',[0.95 0.95 0.1], 'EdgeColor','none', 'FaceAlpha',0.45);
            end
        end
        axis equal vis3d; camlight headlight; lighting gouraud;
        grid on; box on;
        xlabel('X_RAS (mm)'); ylabel('Y_RAS (mm)'); zlabel('Z_RAS (mm)');
        title('metal & ROI, estimated electrode location');
        view(135,20);
        filename = [subj_id, '_identifying_electrodes_on_lead_', traj_id];
        full_path = fullfile(report_dir, [filename, '.pdf']);
        saveas(gcf, full_path, 'pdf');

        % slice plot
        fov_mm=110; res_mm=0.2;
        for k=1:size(contact_centers, 1)
            [I2,Xmm,Ymm] = reslice_plane_RAS(contact_centers(k,:), contact_tangents(k,:), Avox2ras0, ctVol, fov_mm, res_mm);
            figure('Color','w'); imagesc(Xmm(1,:), Ymm(:,1), I2); axis image ij; colormap gray; colorbar;
            hold on; th=linspace(0,2*pi,128);
            plot(rad_mm*cos(th), rad_mm*sin(th), 'y-', 'LineWidth',2);
            plot(0,0,'r+','MarkerSize',10,'LineWidth',1.5);
            title(sprintf('Contact %d (edge-based offsets) — brightest overlap?', k));
            xlabel('v1 (mm)'); ylabel('v2 (mm)');
            filename = [subj_id, 'identifying_electrodes_on_lead_', traj_id, '_CT_electrode_', sprintf('%02d', k), '1_being_distal'];
            full_path = fullfile(report_dir, [filename, '.pdf']);
            saveas(gcf, full_path, 'pdf');
        end

        % line: out.centerline_pts

        % lateral adjustment - needed?
        % --- lateral refine each contact (small in-plane offsets) ---
        % lat_win = 0.5;             % mm (radius of search)
        % nLat   = 7;                % grid resolution per axis
        % lat_vals = linspace(-lat_win, +lat_win, nLat);
        % 
        % for k = 1:size(contact_centers,1)
        %     Ck = contact_centers(k,:); uk = contact_tangents(k,:);
        %     [v1,v2,~] = axis_frame(uk);
        % 
        %     best_mu = -inf; best_C = Ck;
        %     for dx = lat_vals
        %         for dy = lat_vals
        %             Ccand = Ck + dx*v1 + dy*v2;
        %             mu = mean_intensity_cylinder_RAS(Ccand, uk, rad_mm, len_mm, Avox2ras, double(ctVol));
        %             if mu > best_mu
        %                 best_mu = mu; best_C = Ccand;
        %             end
        %         end
        %     end
        %     contact_centers(k,:)  = best_C;
        %     contact_intensity(k)  = best_mu;   % updated score
        %     % (optionally recompute uk via tan_at_s if you can map best_C back to s)
        % end
        
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
        table_filename = [imaging_processing_temp_folder, ...
            '/contact_data_table.mat'];
        save(table_filename, 'contact_tbl', '-v7.3');
        contact_data_struct = table2struct(contact_tbl);
        table_filename = [imaging_processing_temp_folder, ...
            '/contact_data_table_struct.mat'];
        save(table_filename, 'contact_data_struct', '-v7.3');
        imaging_processing_info_struct.shank_model_all = shank_model_all;
        save(updated_imaging_processing_info_struct_save_path, 'imaging_processing_info_struct', '-v7.3');

    end
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
    fid = fopen([automated_aligned_to_brightest_voxel_macro_contacts_folder, '/', file_name], 'w');
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
    fid = fopen([automated_shank_center_line_folder, '/', file_name], 'w');
    fprintf(fid, '%.6f %.6f %.6f\n', C.'); 
    fprintf(fid, 'info\nnumpoints %d\nuseRealRAS 1\n', size(C,1)); fclose(fid);
end



%% get micro
for i_research_shank = 1:length(research_shank_id)
    current_research_shank = research_shank_id{i_research_shank};
    filtered_table = contact_tbl(strcmp(contact_tbl.ShankID, ...
        current_research_shank), :);
    % get coordinates along the shank, fit a line, estimate the micro
    % contact
    x_coords = filtered_table.X(:);
    y_coords = filtered_table.Y(:);
    z_coords = filtered_table.Z(:);
    A = [x_coords y_coords z_coords];
    meanPoint = mean(A, 1);
    % Perform singular value decomposition
    % explain here: https://ltetrel.github.io/data-science/2018/06/08/line_svd.html
    [U,S,V] = svd(A - meanPoint);
    % it's good to know that Ax+By+Cz=D defines a plane. 
    % 3d line: (x,y,z) = (x0,y0,z0)+t(a,b,c)
    % Extract directional vector
    direction = V(:,1);
    first_macro = filtered_table((filtered_table.ContactID == 1), :);
    first_macro_x_ras = first_macro.X;
    first_macro_y_ras = first_macro.Y;
    first_macro_z_ras = first_macro.Z;
    % we use minus because we want contact 8 to contact 1 direction
    estimated_micro_x_ras = first_macro_x_ras - micro_contact_from_first_macro * direction(1);
    estimated_micro_y_ras = first_macro_y_ras - micro_contact_from_first_macro * direction(2);
    estimated_micro_z_ras = first_macro_z_ras - micro_contact_from_first_macro * direction(3);

    estimated_micro_voxel_coords = imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras1 \ [estimated_micro_x_ras;...
    estimated_micro_y_ras; estimated_micro_z_ras; 1];
    first_macro_voxel_coords = imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras1 \ [first_macro_x_ras;...
    first_macro_y_ras; first_macro_z_ras; 1];

    % get a subvolume for micro search
    col_voxel = estimated_micro_voxel_coords(1);
    row_voxel = estimated_micro_voxel_coords(2);
    slice_voxel = estimated_micro_voxel_coords(3);
    % start from 2*2*2 volume, gradually increase volume, get intensity
    % threshold for each search, we will pick the search with highest
    % intensity threshold
    max_search_volume_dim = round(micro_contact_from_first_macro / ...
    max([imaging_processing_info_struct.ct_volume_registered_to_mri.xsize, ...
        imaging_processing_info_struct.ct_volume_registered_to_mri.ysize, imaging_processing_info_struct.ct_volume_registered_to_mri.zsize]));
    intensity_for_each_search = zeros(max_search_volume_dim, 1);
    centroid_coordinates_ras_space_for_each_search = {};
    for search_dim = 1:max_search_volume_dim
        % if sub_vol_voxel_semi_range == 2, then it is 5*5*5 volume search
    % 2*2*2 search Define the neighborhood around the original coordinates
        col_voxel_range = max(1, round(col_voxel - search_dim)):...
            min(size(imaging_processing_info_struct.ct_volume_registered_to_mri.vol, 2), round(col_voxel + search_dim));
        row_voxel_range = max(1, round(row_voxel - search_dim)):...
            min(size(imaging_processing_info_struct.ct_volume_registered_to_mri.vol, 1), round(row_voxel + search_dim));
        slice_voxel_range = max(1, round(slice_voxel - search_dim)):...
            min(size(imaging_processing_info_struct.ct_volume_registered_to_mri.vol, 3), round(slice_voxel + search_dim));
            % Extract the subvolume within the defined neighborhood
        subvolume_manually_identified_contact = imaging_processing_info_struct.ct_volume_registered_to_mri.vol(row_voxel_range, ...
            col_voxel_range, slice_voxel_range);
    
        [col_indices_subvolume_manually_identified_contact, ...
            row_indices_subvolume_manually_identified_contact, ...
            slice_indices_subvolume_manually_identified_contact] = ...
            meshgrid(col_voxel_range, ...
            row_voxel_range, slice_voxel_range);

        [f,xi] = ksdensity(subvolume_manually_identified_contact(:));

        [noisy_voxel_peaks, noisy_voxel_locs] = findpeaks(f);
        last_peak_loc = noisy_voxel_locs(end);
        intensity_threshold = xi(last_peak_loc); % You can set your own threshold
        found_connected_volume = false;

        while ~found_connected_volume
            binary_volume = subvolume_manually_identified_contact > intensity_threshold;
            conn_comp = bwconncomp(binary_volume);
            if (sum(binary_volume(:)) <= minimum_voxel_for_micro_contact) || ...
                (conn_comp.NumObjects > 1)
                intensity_threshold = intensity_threshold - intensity_search_resolution;
            else
                found_connected_volume = true;
                intensity_for_each_search(search_dim) = intensity_threshold;
            end
        end
        % get centroid of the identified volume
        labeled_clusters = bwlabeln(binary_volume);
        cluster_props = regionprops3(labeled_clusters, 'Centroid', 'Volume');
        % col, row, slice, ready to be multiplied by vox2ras1
        centroid_coordinates_voxel_space = [
        cluster_props.Centroid(1) + min(col_indices_subvolume_manually_identified_contact(:)) - 1, ...
        cluster_props.Centroid(2) + min(row_indices_subvolume_manually_identified_contact(:)) - 1, ...
        cluster_props.Centroid(3) + min(slice_indices_subvolume_manually_identified_contact(:)) - 1];
    
        centroid_coordinates_ras_space = imaging_processing_info_struct.ct_volume_registered_to_mri.vox2ras1 * ...
            [centroid_coordinates_voxel_space 1]';
        centroid_coordinates_ras_space = centroid_coordinates_ras_space(1:3);
        centroid_coordinates_ras_space_for_each_search{search_dim} = centroid_coordinates_ras_space;
    end
    % the idea is that, if we pick up scatters from macro contact, we will
    % need a lower intensity threshold to make it connected
    max_intensity_threshold_id = find(intensity_for_each_search == max(intensity_for_each_search));
    micro_coordinates_ras_space = centroid_coordinates_ras_space_for_each_search( ...
        max_intensity_threshold_id);
    micro_coordinates_ras_space = micro_coordinates_ras_space{1};
    % export microelectrode .dat file
    file_name = sprintf('%s_micro.dat', current_research_shank);
    fid = fopen([aligned_to_brightest_voxel_macro_contacts_folder, '/', file_name], 'w');
    fprintf(fid, '\n');
    fprintf(fid, [num2str(micro_coordinates_ras_space') '\n']);
    fprintf(fid, 'info\n');
    fprintf(fid, 'numpoints 1\n');
    fprintf(fid, 'useRealRAS 1\n');
    fclose(fid);
    % 99 representing micro
    new_row = table({current_research_shank}, micro_contact_id, nan, nan, nan, ...
        micro_coordinates_ras_space(1), micro_coordinates_ras_space(2), ...
        micro_coordinates_ras_space(3), 'VariableNames', ...
        contact_tbl.Properties.VariableNames);
    contact_tbl = [contact_tbl; new_row];
end
