function model = dixi_model_from_row(n_contact_str, n_contact_numeric)
% Parse DIXI depth electrode model and generate offsets from distal (mm). updated 2025, setting: BJH, st. Louis
% Inputs:
%   n_contact_str     e.g., 'D08-15AM', 'D08-15BM', 'D08-15CM', 'D08-18CM', or '15 (3x5)'
%   n_contact_numeric numeric fallback if string missing (e.g., 15, 18)
% Output struct:
%   .type            'AM'|'BM'|'CM'
%   .n_contact       scalar
%   .groups          [G, K] if grouped; [] for AM
%   .intra_pitch_mm  within-group pitch (mm)
%   .inter_gap_mm    between-group gap (mm) (AM: 3.5)
%   .diameter_mm     contact diameter (default 0.8)
%   .contact_len_mm  contact length (default 2.0)
%   .offsets_mm      1×N offsets from distal tip (mm), ascending

    % Defaults
    model = struct('type','AM','n_contact',[], 'groups',[], ...
                   'intra_pitch_mm',3.5, 'inter_gap_mm',3.5, ...
                   'diameter_mm',0.8, 'contact_len_mm',2.0, ...
                   'offsets_mm',[]);
    % Normalize inputs
    str = '';
    if ~isempty(n_contact_str) && ischar(n_contact_str) || (isstring(n_contact_str) && strlength(n_contact_str)>0)
        str = char(n_contact_str);
    end

    % Heuristics to detect type
    isAM = contains(str,'AM','IgnoreCase',true);
    isBM = contains(str,'BM','IgnoreCase',true);
    isCM = contains(str,'CM','IgnoreCase',true);

    % Extract N and group pattern if present
    N = [];
    tokN = regexp(str, '(?<!\d)(\d+)(?!\d)', 'tokens', 'once'); % first integer
    if ~isempty(tokN), N = str2double(tokN{1}); end

    grp = []; % [GxK]
    tokGXK = regexp(str, '(\d+)\s*\(\s*(\d+)\s*[xX]\s*(\d+)\s*\)', 'tokens', 'once'); % e.g., '15 (3x5)'
    if ~isempty(tokGXK)
        % The middle number is total contacts, ignore; we only need GxK
        G = str2double(tokGXK{2});
        K = str2double(tokGXK{3});
        grp = [G, K];
    else
        % Infer standard patterns for 15/18 when BM/CM
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
            error('Cannot determine n_contact.');
        end
    end

    % Decide type
    if isBM, model.type='BM'; end
    if isCM, model.type='CM'; end
    if ~(isBM||isCM), model.type='AM'; end

    % Geometry by type
    switch model.type
        case 'AM'
            model.n_contact = N;
            model.intra_pitch_mm = 3.5;
            model.inter_gap_mm   = 3.5; % same as intra (uniform)
            model.groups = [];
            model.offsets_mm = (0:N-1) * model.intra_pitch_mm;

        case {'BM','CM'}
            % Required group pattern
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

            % Offsets: group blocks spaced by inter_gap; within group spaced by intra_pitch
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

        otherwise
            error('Unknown model type.');
    end
end

function y = iff(cond, a, b)
    if cond, y=a; else, y=b; end
end