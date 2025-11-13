function model = adtech_model_from_row(n_contact_str, n_contact_numeric)
%ADTECH_MODEL_FROM_ROW Generate a Behnke-Fried (BF) electrode model for Ad-Tech.
%
% This function returns a geometry model for Ad-Tech BF electrodes, which consist of
% 8 macro contacts arranged linearly with 5 mm center-to-center spacing.
%
% INPUTS:
%   n_contact_str      : (optional) A string or character description of the contact model (e.g., 'BF')
%                        This input is currently ignored but included for interface compatibility.
%   n_contact_numeric  : (optional) Numeric count of macro contacts (e.g., 8). If provided and not 8,
%                        a warning is issued.
%
% OUTPUT:
%   model : Struct containing the following fields:
%       - type               : Electrode type ('BF')
%       - n_contact          : Number of macro contacts (expected to be 8)
%       - groups             : [] (no grouping for BF)
%       - intra_pitch_mm     : 5.0 (spacing between adjacent contacts)
%       - inter_gap_mm       : 5.0 (same as intra_pitch, as BF has no groupings)
%       - diameter_mm        : 1.28 (contact diameter in mm)
%       - contact_len_mm     : 1.57 (contact length in mm)
%       - offsets_mm         : [0 5 10 15 20 25 30 35] (contact center offsets from distal tip, mm)
%       - inactive_idx       : [] (all contacts are active)
%       - active_idx         : 1:8 (indices of active contacts)
%       - active_offsets_mm  : same as offsets_mm (since all contacts are active)
%
% The output is compatible with the structure returned by `dixi_model_from_row`.

    % Define fixed geometry for Ad-Tech BF electrodes
    contact_count = 8;
    pitch_mm = 5.0;
    offsets = 0:pitch_mm:(pitch_mm * (contact_count - 1));  % [0 5 10 ... 35]

    model = struct( ...
        'type', 'BF', ...
        'n_contact', contact_count, ...
        'groups', [], ...
        'intra_pitch_mm', pitch_mm, ...
        'inter_gap_mm', pitch_mm, ...  % No grouping → same as intra
        'diameter_mm', 1.28, ...
        'contact_len_mm', 1.57, ...
        'offsets_mm', offsets, ...
        'inactive_idx', [], ...
        'active_idx', 1:contact_count, ...
        'active_offsets_mm', offsets ...
    );

    % Validate contact count if provided
    if nargin >= 2 && ~isempty(n_contact_numeric)
        if n_contact_numeric ~= contact_count
            warning('Ad-Tech BF electrodes are expected to have 8 contacts, but got %d.', n_contact_numeric);
        end
    end

    % Assertion: Must have exactly 8 active contacts
    assert(numel(model.active_offsets_mm) == 8, ...
        'adtech_model_from_row:InvalidActiveCount', ...
        'Expected 8 active offsets for BF electrode, got %d.', numel(model.active_offsets_mm));
end