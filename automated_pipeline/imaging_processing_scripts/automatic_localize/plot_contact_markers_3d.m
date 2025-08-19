function plot_contact_markers_3d(C, r, rgb, labelStr)
% Plot small spheres at rows of C (Nx3). Labels only the first for legend.
% C: [x y z], r: radius, rgb: 1x3 color, labelStr: e.g., 'orig' or 'rob'

    if isempty(C), return; end
    if nargin < 4, labelStr = ''; end

    try
        [sx,sy,sz] = sphere(12);
        h = gobjects(size(C,1),1);
        for i = 1:size(C,1)
            h(i) = surf(r*sx + C(i,1), r*sy + C(i,2), r*sz + C(i,3), ...
                'EdgeColor','none','FaceColor',rgb,'FaceAlpha',0.9);
        end
        % Legend: name only the first, hide the rest
        if ~isempty(labelStr)
            set(h(1), 'DisplayName', sprintf('centers (%s)', labelStr), 'HandleVisibility','on');
        else
            set(h(1), 'HandleVisibility','on');
        end
        if numel(h) > 1, set(h(2:end), 'HandleVisibility','off'); end
    catch
        % Fallback: points instead of spheres
        h = scatter3(C(:,1), C(:,2), C(:,3), 36, 'MarkerFaceColor',rgb, 'MarkerEdgeColor','k');
        if ~isempty(labelStr)
            set(h, 'DisplayName', sprintf('centers (%s)', labelStr), 'HandleVisibility','on');
        else
            set(h, 'HandleVisibility','on');
        end
    end
end