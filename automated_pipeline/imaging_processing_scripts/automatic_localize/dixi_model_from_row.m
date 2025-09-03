function model = dixi_model_from_row(n_contact_str, n_contact_numeric)
% Parse DIXI depth electrode model and generate offsets from distal (mm).
% Inputs:
%   n_contact_str     e.g., 'D08-15AM', 'D08-15BM', 'D08-15CM', 'D08-18CM',
%                      '15 (3x5)', 'MME-11', 'MME'
%   n_contact_numeric numeric fallback if string missing (e.g., 15, 18)
% Output struct:
%   .type               'AM'|'BM'|'CM'|'MME'
%   .n_contact          scalar (total, incl. inactive for MME)
%   .groups             [G, K] if grouped; [] otherwise
%   .intra_pitch_mm     within-group pitch (mm) or uniform pitch
%   .inter_gap_mm       between-group gap (mm) (AM: 3.5)
%   .diameter_mm        contact diameter (default 0.8)
%   .contact_len_mm     contact length (default 2.0)
%   .offsets_mm         1×N offsets from distal tip (mm), ascending
%   .inactive_idx       indices of contacts that are not connected (MME only)
%   .active_idx         indices of contacts that are connected (MME only)
%   .active_offsets_mm  offsets for active contacts (MME only)

    % Defaults, % DIXI models all have 0.8mm diameter
    model = struct('type','AM','n_contact',[], 'groups',[], ...
                   'intra_pitch_mm',3.5, 'inter_gap_mm',3.5, ...
                   'diameter_mm',0.8, 'contact_len_mm',2.0, ...
                   'offsets_mm',[], ...
                   'inactive_idx',[], 'active_idx',[], 'active_offsets_mm',[]);
    % Normalize inputs
    str = '';
    if (~isempty(n_contact_str) && (ischar(n_contact_str) || isstring(n_contact_str))) ...
            && strlength(string(n_contact_str))>0
        str = char(n_contact_str);
    end

    % Heuristics to detect type
    isAM  = contains(str,'AM','IgnoreCase',true);
    isBM  = contains(str,'BM','IgnoreCase',true);
    isCM  = contains(str,'CM','IgnoreCase',true);
    isMME = contains(str,'MME','IgnoreCase',true) || ...
            contains(str,'macro micro','IgnoreCase',true);

    % Extract N and group pattern if present
    N = [];
    tokN = regexp(str, '(?<!\d)(\d+)(?!\d)', 'tokens', 'once'); % first integer
    if ~isempty(tokN), N = str2double(tokN{1}); end

    grp = []; % [GxK]
    tokGXK = regexp(str, '(\d+)\s*\(\s*(\d+)\s*[xX]\s*(\d+)\s*\)', 'tokens', 'once'); % e.g., '15 (3x5)'
    if ~isempty(tokGXK)
        G = str2double(tokGXK{2});
        K = str2double(tokGXK{3});
        grp = [G, K];
    else
        if (isBM || isCM) && ~isempty(N)
            if N==15, grp = [3,5]; end
            if N==18, grp = [3,6]; end
        end
    end

    % Fallback N if missing
    if isempty(N)
        if ~isempty(n_contact_numeric) && ~isnan(n_contact_numeric)
            N = n_contact_numeric;
        else
            % For MME default to 11 if nothing specified
            if isMME, N = 11; else, error('Cannot determine n_contact.'); end
        end
    end

    % Decide type
    if isBM,  model.type='BM'; end
    if isCM,  model.type='CM'; end
    if isMME, model.type='MME'; end
    if ~(isBM||isCM||isMME), model.type='AM'; end

    % Geometry by type
    switch model.type
        case 'AM'
            model.n_contact = N;
            model.intra_pitch_mm = 3.5;
            model.inter_gap_mm   = 3.5; % same as intra (uniform)
            model.groups = [];
            model.offsets_mm = (0:N-1) * model.intra_pitch_mm;

        case {'BM','CM'}
            if isempty(grp)
                error('Grouped model %s but no GxK parsed (e.g., 3x5).', model.type);
            end
            G = grp(1); K = grp(2);
            if G*K ~= N
                warning('Group pattern %dx%d != n_contact %d. Adjusting N to %d.', G, K, N, G*K);
                N = G*K;
            end
            model.n_contact = N;
            model.groups = [G, K];
            model.intra_pitch_mm = 1.5;
            model.inter_gap_mm   = iff(strcmpi(model.type,'BM'), 7.0, 11.0);

            offsets = zeros(1, N);
            idx = 1;
            for g = 1:G
                base_g = (g-1) * ( (K-1)*model.intra_pitch_mm + model.inter_gap_mm );
                for k = 1:K
                    offsets(idx) = base_g + (k-1)*model.intra_pitch_mm;
                    idx = idx + 1;
                end
            end
            model.offsets_mm = offsets;

        case 'MME'
            % Your MME definition:
            % - 2 mm contact length
            % - 4 mm inter-contact spacing (uniform pitch of 2 mm)
            % - 11 contacts total
            % - middle 5 contacts are NOT connected to amplifier (inactive), but still visible on CT
             % inter_gap_mm:
	        %Definition: extra spacing between groups (for BM/CM types)
            model.n_contact = N;
            model.contact_len_mm = 2.0;
            model.intra_pitch_mm = 4.0;   % center-to-center pitch
            model.inter_gap_mm   = 2.0;   % uniform
            model.groups = [];            % single row
            total_n_contact_in_CT = N + 5;
            model.offsets_mm = (0:total_n_contact_in_CT-1) * model.intra_pitch_mm;

            % Define middle-5 inactive for any odd N>=5 (centered set)
            
            if total_n_contact_in_CT >= 5
                mid  = ceil(total_n_contact_in_CT/2);       % center index (for N=11 -> 6)
                half = 2;               % (5-1)/2
                inactive = (mid-half):(mid+half); % 5 indices (for 11 -> 4:8)
            else
                inactive = [];
            end
            model.inactive_idx = inactive(:).';
            model.active_idx   = setdiff(1:total_n_contact_in_CT, model.inactive_idx);
            model.active_offsets_mm  = model.offsets_mm(model.active_idx);

        otherwise
            error('Unknown model type.');
    end
end

function y = iff(cond, a, b)
    if cond, y=a; else, y=b; end
end