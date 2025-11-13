function model = pmt_model_from_row(n_contact_str, n_contact_numeric)
%PMT_MODEL_FROM_ROW Generate a PMT depth electrode model structure.
%
% Handles standard PMT depth electrodes with variable contact counts
% and spacing (e.g., 3.5 mm default, 3.97 mm for '16-092', 4.43 mm for '16-093').
% Also includes model '16-100' with two 8-contact groups separated by 35 mm.
%
% INPUTS:
%   n_contact_str     : (optional) A string like '8', '16-093', '16-100'
%   n_contact_numeric : (optional) A number like 8 or 16
%
% OUTPUT:
%   model : Struct with standard fields:
%       - type               : Electrode type (e.g., 'PMT', '16-093')
%       - n_contact          : Number of macro contacts
%       - groups             : Cell array of contact indices per group
%       - intra_pitch_mm     : Spacing between contacts (in mm)
%       - inter_gap_mm       : Gap between groups (if applicable)
%       - diameter_mm        : 0.8 mm (typical PMT diameter)
%       - contact_len_mm     : 2 mm (PMT macro contact length)
%       - offsets_mm         : Contact center offsets from tip
%       - inactive_idx       : [] (all active)
%       - active_idx         : 1:N
%       - active_offsets_mm  : same as offsets_mm

    % Defaults
    contact_diameter_mm = 0.8;
    contact_length_mm   = 2;
    default_spacing     = 3.5;
    default_gap         = default_spacing;

    % Parse model type and spacing
    if ischar(n_contact_str) || isstring(n_contact_str)
        model_type = char(n_contact_str);
        switch model_type
            case '16-092'
                pitch_mm = 3.97;
                inter_gap_mm = pitch_mm;
                grouped = false;

            case '16-093'
                pitch_mm = 4.43;
                inter_gap_mm = pitch_mm;
                grouped = false;

            case '16-100'
                pitch_mm = 3.5;
                inter_gap_mm = 36; % between two groups
                grouped = true;

            otherwise
                pitch_mm = default_spacing;
                inter_gap_mm = pitch_mm;
                grouped = false;
                model_type = 'PMT';
        end
    else
        pitch_mm = default_spacing;
        inter_gap_mm = default_gap;
        grouped = false;
        model_type = 'PMT';
    end

    % Infer number of contacts
    if nargin < 2 || isempty(n_contact_numeric)
        error('n_contact_numeric must be provided.');
    end
    n_contact = n_contact_numeric;

    % Compute offsets
    if grouped && strcmp(model_type, '16-100')
        % Two groups of 8 contacts
        group1 = 0:pitch_mm:(pitch_mm * 7);
        group2_start = group1(end) + pitch_mm + inter_gap_mm; % add the 38.5 mm gap
        group2 = group2_start + (0:pitch_mm:(pitch_mm * 7));

        offsets = [group1, group2];
        groups = {1:8, 9:16};
    else
        offsets = 0:pitch_mm:(pitch_mm * (n_contact - 1));
        groups = [];
    end

    % Build model struct
    model = struct( ...
        'type', model_type, ...
        'n_contact', n_contact, ...
        'groups', {groups}, ...
        'intra_pitch_mm', pitch_mm, ...
        'inter_gap_mm', inter_gap_mm, ...
        'diameter_mm', contact_diameter_mm, ...
        'contact_len_mm', contact_length_mm, ...
        'offsets_mm', offsets, ...
        'inactive_idx', [], ...
        'active_idx', 1:n_contact, ...
        'active_offsets_mm', offsets ...
    );

    % Sanity check
    if length(offsets) ~= n_contact
        warning('Mismatch between expected contact count (%d) and offset vector length (%d).', ...
            n_contact, length(offsets));
    end
end