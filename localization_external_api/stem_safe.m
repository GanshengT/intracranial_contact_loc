function h = stem_safe(x, y, varargin)
% STEM_SAFE  Minimal wrapper around stem() that ignores NaN/Inf in y
% Usage:
%   stem_safe(x, y, 'Color',[...], 'Marker','o', 'DisplayName','text', 'LineWidth',1.2, 'MarkerSize',36)

    % defaults
    p = inputParser;  p.CaseSensitive = false;
    addParameter(p,'Color',            [0 0 0]);
    addParameter(p,'Marker',           'o');      % stem supports 'Marker'
    addParameter(p,'DisplayName',      '');
    addParameter(p,'LineWidth',        1.0);
    addParameter(p,'MarkerSize',       10);       % applies to stem markers
    parse(p, varargin{:});
    P = p.Results;

    x = x(:); y = y(:);
    good = ~isnan(y) & ~isinf(y);
    if ~any(good)
        h = gobjects(1);
        return;
    end

    % Make the stem and then set properties (don’t pass PV pairs directly)
    h = stem(x(good), y(good));  % no 'filled' here
    set(h, 'Color', P.Color, 'Marker', P.Marker, 'LineWidth', P.LineWidth);
    % Some MATLAB versions use 'MarkerSize' on the Stem object; apply if present:
    if isprop(h,'MarkerSize'), set(h,'MarkerSize',P.MarkerSize); end
    if ~isempty(P.DisplayName), set(h, 'DisplayName', P.DisplayName); end
end