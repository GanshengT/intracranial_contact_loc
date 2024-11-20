%% Description
% This script localize the micro contact, adjust manually localized macro
% contact
% This script outputs the contact coordinates and label
% Figures showing the contact location will be generated
% 2024 version
% contact: Gansheng Tan g.tan@wustl.edu
clear;
close all;
%% definition
% we load all the info
addpath('/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/projects/BLEAS/theta_saccade/script');
addpath(genpath(['/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/' ...
    'projects/BLEAS/theta_saccade/script/imaging_processing_scripts']));
data_path = '/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/projects/BLEAS/theta_saccade/data';
% subj info
subj_id = 'BJH025';
% research shank name
% this is subject-specific, be sure that you change that. This should match
% manually identified labels
research_shank_id = {'HL_MTG_AHC_BEFR', 'HR_MTG_AHC_BEFR'};
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
% imaging_processing_temp_folder = [data_path, '/', subj_id, '/imaging_process'];
% imaging_processing_info_struct_save_path = [imaging_processing_temp_folder, ...
%     '/imaging_processing_info_struct.mat'];
imaging_processing_info_struct_save_path = 'imaging_processing_info_struct.mat';
imaging_processing_info_struct = load(imaging_processing_info_struct_save_path);
% if you do not have manually identified macro contact ready yet
% uncomment below for manual edition in freeview
% !! important, from proximal to distal, id number goes from 1 to 16 !!
% OR_sheet_path = [data_path, '/' subj_id, '/notes/SEEG-OR-', subj_id, '.xlsx'];
% sheetname = 'Surgeon Date';
% col_name_electrode_number = 'Est_Electrode';
% electrode_info = ...
% extract_electrode_info_fromOR_sheet(OR_sheet_path, sheetname, col_name_electrode_number);
% manual_macro_contact_identification(['/Users/ganshengtan/Library/' ...
%     'CloudStorage/Box-Box/Washu/projects/BLEAS/theta_saccade/data/BJH025/imaging_process'], ...
%     'SOFT_TISSUE_1_0_registered_to_mri', 'SAG_T1_MPRAGE_converted_from_orig_mgz', electrode_info);

% load manually identified macro contact
% manually_identified_macro_contacts_path = [imaging_processing_temp_folder, ...
%     '/manually_identified_macro_contacts/'];
manually_identified_macro_contacts_path = 'manually_identified_macro_contacts/';
% save macro contacts' location that are aligned with the brightest voxel
% aligned_to_brightest_voxel_macro_contacts_folder = [imaging_processing_temp_folder, ...
%     '/aligned_to_brightest_voxel_macro_contacts'];
aligned_to_brightest_voxel_macro_contacts_folder = 'aligned_to_brightest_voxel_macro_contacts';
if ~exist(aligned_to_brightest_voxel_macro_contacts_folder, 'dir')
    mkdir(aligned_to_brightest_voxel_macro_contacts_folder);
end
% save contacts' location visualization - see tutorial in a separate scripts
% report_contacts_related_surf_folder = [imaging_processing_temp_folder, ...
%     '/report_contacts_related_surf_folder'];
% if ~exist(report_contacts_related_surf_folder, 'dir')
%     mkdir(report_contacts_related_surf_folder);
% end

% region of interest
% we will focus on aparc - Desikan/Killiany parcellation

% region_of_interest_to_add = {};
% for i_region_of_interest = 1:length(region_of_interest)
%     current_region = region_of_interest{i_region_of_interest};
%     if contains(current_region, 'Left')
%         new_region = strrep(current_region, 'Left', 'Right');
%         region_of_interest_to_add = [region_of_interest_to_add, new_region];
%     end
% 
%     if contains(current_region, 'lh')
%         new_region = strrep(current_region, 'lh', 'rh');
%         region_of_interest_to_add = [region_of_interest_to_add, new_region];
%     end
% end
% 
% region_of_interest = [region_of_interest, region_of_interest_to_add];
% region_of_interest_to_add = {};
% for i_region_of_interest = 1:length(region_of_interest_hippo_subsegment)
%     current_region = region_of_interest_hippo_subsegment{i_region_of_interest};
%     if contains(current_region, 'head')
%         new_region = strrep(current_region, 'head', 'body');
%         region_of_interest_to_add = [region_of_interest_to_add, new_region];
%     end
% end
% region_of_interest = [region_of_interest, region_of_interest_to_add];

% create struct that stores contact coordinates
contact_info = struct();
% contact_info stores shank_id, contact_id, contact_type (macro or micro),
% coordinates, and label.

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


%% adjust macro contact
ct_volume_registered_to_mri = imaging_processing_info_struct.ct_volume_registered_to_mri;
% voxel_volume
voxel_volume = ct_volume_registered_to_mri.xsize * ...
ct_volume_registered_to_mri.ysize * ct_volume_registered_to_mri.zsize;
minimum_voxel_for_contact = ceil(macro_contact_volume_mm_cube / voxel_volume);

% we will first get RAS coordinates corresponding to ct_volume_registered_to_mri.vol
ct_volume_registered_to_mri_ras_coords = get_ras_coordinates_from_vox(ct_volume_registered_to_mri);
% in freeview, it displays voxel coordinate as [col, row, slice], each
% starting from 0 0 0
% ras_coordinates of [100, 1, 100] voxel space in freeview correspond to 
% ras_coordinates of ct_volume_registered_to_mri_ras_coords(2, 101, 101)
% ct_volume_registered_to_mri.vol(2, 101, 101) = freeview_ct_vol(100,1,100)
% so the voxel value of ct_volume_registered_to_mri_ras_coords(2, 101, 101)
% = ct_volume_registered_to_mri.vol(2, 101, 101)
sub_vol_voxel_range = ceil(sub_vol_range_mm / min(ct_volume_registered_to_mri.volres));
sub_vol_voxel_semi_range = ceil(sub_vol_voxel_range / 2);

% Loop through each row in contact_data_table
for i_row = 1:size(contact_data_table, 1)
    % Extract the original coordinates from the table
    manually_identified_contact_x_ras = contact_data_table.X_manual(i_row);
    manually_identified_contact_y_ras = contact_data_table.Y_manual(i_row);
    manually_identified_contact_z_ras = contact_data_table.Z_manual(i_row);

    voxel_coords = ct_volume_registered_to_mri.vox2ras1 \ [manually_identified_contact_x_ras;...
        manually_identified_contact_y_ras; manually_identified_contact_z_ras; 1];

    col_voxel = round(voxel_coords(1));
    row_voxel = round(voxel_coords(2));
    slice_voxel = round(voxel_coords(3));
    
    % Define the neighborhood around the original coordinates
    col_voxel_range = max(1, col_voxel - sub_vol_voxel_semi_range):...
        min(size(ct_volume_registered_to_mri.vol, 2), col_voxel + sub_vol_voxel_semi_range);
    row_voxel_range = max(1, row_voxel - sub_vol_voxel_semi_range):...
        min(size(ct_volume_registered_to_mri.vol, 1), row_voxel + sub_vol_voxel_semi_range);
    slice_voxel_range = max(1, slice_voxel - sub_vol_voxel_semi_range):...
        min(size(ct_volume_registered_to_mri.vol, 3), slice_voxel + sub_vol_voxel_semi_range);

    % Extract the subvolume within the defined neighborhood
    subvolume_manually_identified_contact = ct_volume_registered_to_mri.vol(row_voxel_range, ...
        col_voxel_range, slice_voxel_range);

    [col_indices_subvolume_manually_identified_contact, ...
        row_indices_subvolume_manually_identified_contact, ...
        slice_indices_subvolume_manually_identified_contact] = ...
        meshgrid(col_voxel_range, ...
        row_voxel_range, slice_voxel_range);

    % get cdf of voxel value and make cut might not be good idea, unstable
    % method
    % DBSCAN is not good neither because the intensity of contact voxel varies
    % here we develop an intensity-based method
    % estimate kernel density
    [f,xi] = ksdensity(subvolume_manually_identified_contact(:));
    % find the first peak (noisy voxels, not contacts), and start our
    % search from that value
    [noisy_voxel_peaks, noisy_voxel_locs] = findpeaks(f);
    last_peak_loc = noisy_voxel_locs(end);
    intensity_threshold = xi(last_peak_loc); % You can set your own threshold
    found_connected_volume = false;

    while ~found_connected_volume
        binary_volume = subvolume_manually_identified_contact > intensity_threshold;
        conn_comp = bwconncomp(binary_volume);
        if (sum(binary_volume(:)) <= minimum_voxel_for_contact) || ...
                (conn_comp.NumObjects > 1)
            intensity_threshold = intensity_threshold - intensity_search_resolution;
        else
            found_connected_volume = true;
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

    centroid_coordinates_ras_space = ct_volume_registered_to_mri.vox2ras1 * ...
        [centroid_coordinates_voxel_space 1]';
    centroid_coordinates_ras_space = centroid_coordinates_ras_space(1:3);
    % add columns to the table
    contact_data_table.X_aligned_to_brightest_voxel(i_row) = centroid_coordinates_ras_space(1);
    contact_data_table.Y_aligned_to_brightest_voxel(i_row) = centroid_coordinates_ras_space(2);
    contact_data_table.Z_aligned_to_brightest_voxel(i_row) = centroid_coordinates_ras_space(3);
end

% export the realigned electrode location and check in freeview
unique_shankIDs = unique(contact_data_table.ShankID);
for i_shankID = 1:length(unique_shankIDs)
    shankID = unique_shankIDs{i_shankID};
    filtered_table = contact_data_table(strcmp(contact_data_table.ShankID, shankID), :);
    X = filtered_table.X_aligned_to_brightest_voxel;
    Y = filtered_table.Y_aligned_to_brightest_voxel;
    Z = filtered_table.Z_aligned_to_brightest_voxel;
    file_name = sprintf('%s.dat', shankID);

    fid = fopen([aligned_to_brightest_voxel_macro_contacts_folder, '/', file_name], 'w');
    fprintf(fid, '\n');
    for i_contact = 1:height(filtered_table)
        fprintf(fid, [num2str([X(i_contact), Y(i_contact), Z(i_contact)]) '\n']);
    end
    fprintf(fid, 'info\n');
    fprintf(fid, 'numpoints %d\n', height(filtered_table));
    fprintf(fid, 'useRealRAS 1\n');
    fclose(fid);
end

%% macro objectively identified, get micro
% next we want to get the micro electrode location
% The idea is to fit a line, and extend the line and form a search volume
% then we apply to same technique to localize the micro electrodes.
for i_research_shank = 1:length(research_shank_id)
    current_research_shank = research_shank_id{i_research_shank};
    filtered_table = contact_data_table(strcmp(contact_data_table.ShankID, ...
        current_research_shank), :);
    % get coordinates along the shank, fit a line, estimate the micro
    % contact
    x_coords = filtered_table.X_aligned_to_brightest_voxel(:);
    y_coords = filtered_table.Y_aligned_to_brightest_voxel(:);
    z_coords = filtered_table.Z_aligned_to_brightest_voxel(:);
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
    first_macro_x_ras = first_macro.X_aligned_to_brightest_voxel;
    first_macro_y_ras = first_macro.Y_aligned_to_brightest_voxel;
    first_macro_z_ras = first_macro.Z_aligned_to_brightest_voxel;
    % we use minus because we want contact 8 to contact 1 direction
    estimated_micro_x_ras = first_macro_x_ras - micro_contact_from_first_macro * direction(1);
    estimated_micro_y_ras = first_macro_y_ras - micro_contact_from_first_macro * direction(2);
    estimated_micro_z_ras = first_macro_z_ras - micro_contact_from_first_macro * direction(3);

    estimated_micro_voxel_coords = ct_volume_registered_to_mri.vox2ras1 \ [estimated_micro_x_ras;...
    estimated_micro_y_ras; estimated_micro_z_ras; 1];
    first_macro_voxel_coords = ct_volume_registered_to_mri.vox2ras1 \ [first_macro_x_ras;...
    first_macro_y_ras; first_macro_z_ras; 1];

    % get a subvolume for micro search
    col_voxel = estimated_micro_voxel_coords(1);
    row_voxel = estimated_micro_voxel_coords(2);
    slice_voxel = estimated_micro_voxel_coords(3);

    % start from 2*2*2 volume, gradually increase volume, get intensity
    % threshold for each search, we will pick the search with highest
    % intensity threshold
    max_search_volume_dim = micro_contact_from_first_macro / ...
    max([ct_volume_registered_to_mri.xsize, ...
        ct_volume_registered_to_mri.ysize, ct_volume_registered_to_mri.zsize]);
    intensity_for_each_search = zeros(max_search_volume_dim, 1);
    centroid_coordinates_ras_space_for_each_search = {};
    for search_dim = 1:max_search_volume_dim
        % if sub_vol_voxel_semi_range == 2, then it is 5*5*5 volume search
    % 2*2*2 search Define the neighborhood around the original coordinates
        col_voxel_range = max(1, round(col_voxel - search_dim)):...
            min(size(ct_volume_registered_to_mri.vol, 2), round(col_voxel + search_dim));
        row_voxel_range = max(1, round(row_voxel - search_dim)):...
            min(size(ct_volume_registered_to_mri.vol, 1), round(row_voxel + search_dim));
        slice_voxel_range = max(1, round(slice_voxel - search_dim)):...
            min(size(ct_volume_registered_to_mri.vol, 3), round(slice_voxel + search_dim));
            % Extract the subvolume within the defined neighborhood
        subvolume_manually_identified_contact = ct_volume_registered_to_mri.vol(row_voxel_range, ...
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
    
        centroid_coordinates_ras_space = ct_volume_registered_to_mri.vox2ras1 * ...
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
        contact_data_table.Properties.VariableNames);
    contact_data_table = [contact_data_table; new_row];
end

% contact_data_table is what you needed

