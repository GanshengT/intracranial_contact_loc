function [fill_new, dstar_new, ctr_new, t0_new, t1_new, r_enclose_new] = eval_fill_after_add( ...
        G, cand_ijk, vx, U_in, TPadMM, ...
        cand_base, Xras, Avox2ras0, varargin)
% EVAL_FILL_AFTER_ADD (enclosing-tube version)
% After tentatively adding cand_ijk to G:
%   1) recompute centroid
%   2) refit axis via fit_axis_by_fill (same options you use for cur_fill)
%   3) get enclosing tube via enclosing_tube_params
%   4) build tube mask and compute fill = (|G2∩tube|*Vvoxel) / V_tube
%
% Returns fill_new and the fitted params so caller can reuse if desired.

    % parse and forward options to fit_axis_by_fill
    ip = inputParser; ip.CaseSensitive=false;
    addParameter(ip,'AngleDeg',5,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=45);
    addParameter(ip,'RingK',12,@(x)isnumeric(x)&&isscalar(x)&&x>=4&&x<=128);
    addParameter(ip,'IncludeCenter',true,@(x)islogical(x)||ismember(x,[0 1]));
    addParameter(ip,'RGuardMM',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'TGuardMM',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'RobustPct',98,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=100);
    parse(ip,varargin{:});

    % 0) Copy & add candidate (bounds guard)
    G2 = G;
    i = cand_ijk(1); j = cand_ijk(2); k = cand_ijk(3);
    if i<1 || j<1 || k<1 || i>size(G,1) || j>size(G,2) || k>size(G,3)
        % out of bounds -> no change from current; compute consistent fill anyway
        G2 = G;  % (falls through)
    else
        G2(i,j,k) = true;
    end

    % 1) centroid (RAS)
    ctr_new = blob_centroid_ras(G2, Xras);

    % 2) axis fit (volume-fill, same options as cur_fill’s dstar)
    [dstar_new, ~, ~, ~, ~] = fit_axis_by_fill( ...
        U_in, ctr_new, TPadMM, G2, Xras, vx, ...
        'AngleDeg',     ip.Results.AngleDeg, ...
        'RingK',        ip.Results.RingK, ...
        'IncludeCenter',ip.Results.IncludeCenter, ...
        'RGuardMM',     ip.Results.RGuardMM, ...
        'TGuardMM',     ip.Results.TGuardMM, ...
        'RobustPct',    ip.Results.RobustPct);

    % 3) enclosing tube around G2 given dstar_new (consistent with cur_fill)
    [t0_new, t1_new, r_enclose_new, stats_new] = enclosing_tube_params(G2, Xras, ctr_new, dstar_new); 

    % 4) build tube mask on the SAME domain policy as cur_fill
    %    NOTE: use Domain=(cand_base | G2) to match your current code’s behavior.
    Vvoxel      = prod(vx);                                  % mm^3 per voxel
    V_tube_mm3  = max(1e-9, pi * (r_enclose_new^2) * max(0, (t1_new - t0_new)));
    fill_new    = safe_ratio(nnz(G2) * Vvoxel, V_tube_mm3);

    % diagonostic plot
    % figure('Color','w'); hold on;
    % % Blob surface (G) in RAS
    % if any(G2(:))
    %     [m,n,p] = size(G2);
    %     [Xm,Ym,Zm] = meshgrid(1:n,1:m,1:p);
    %     hg = patch(isosurface(Xm,Ym,Zm,G2,0.5));
    %     set(hg,'EdgeColor','none','FaceColor',[1 0.7 0],'FaceAlpha',0.3);
    %     V0 = hg.Vertices - 1; Vras = [V0, ones(size(V0,1),1)] * Avox2ras0.'; 
    %     hg.Vertices = Vras(:,1:3);
    % end
    % 
    % % Planned line (span length based on voxel bbox extent)
    % 
    % 
    % ctr = blob_centroid_ras(G2, Xras); % your existing helper
    % 
    % % [dstar, cur_fill, t0, t1, r_mm] = fit_axis_by_fill( ...
    % %     U, ctr, TPadMM, G2, Xras, vx, ...
    % %     'AngleDeg', AngleDeg, 'RingK', 12, 'IncludeCenter', true, ...
    % %     'RGuardMM', 0.2, 'TGuardMM', 0.5, 'RobustPct', 98);
    % % 
    % % [t0, t1, r_mm, stats] = enclosing_tube_params(G, Xras, ctr, dstar);
    % 
    % scatter3(ctr(1),ctr(2),ctr(3),80,[0.85 0.33 0.10],'filled','MarkerEdgeColor','k');
    % 
    % % longitudinal sampling (denser if longer tube)
    % L  = max(t1_new - t0_new, 1e-3);
    % nTau = max(24, ceil(L/3));     % ~3 mm spacing along axis
    % nTh  = 48;                      % circular resolution
    % tmp = [1;0;0];
    % if abs(dot(tmp,dstar_new)) > 0.99
    %     tmp = [0;1;0];
    % end
    % n1 = cross(dstar_new, tmp); n1 = n1 / norm(n1);
    % n2 = cross(dstar_new, n1);  n2 = n2 / norm(n2);
    % tau   = linspace(t0_new, t1_new, nTau);
    % theta = linspace(0, 2*pi, nTh);
    % [TH, TAU] = meshgrid(theta, tau);
    % 
    % % parametric tube surface in RAS
    % Xt = ctr(1) + TAU.*dstar_new(1) + r_enclose_new*(cos(TH)*n1(1) + sin(TH)*n2(1));
    % Yt = ctr(2) + TAU.*dstar_new(2) + r_enclose_new*(cos(TH)*n1(2) + sin(TH)*n2(2));
    % Zt = ctr(3) + TAU.*dstar_new(3) + r_enclose_new*(cos(TH)*n1(3) + sin(TH)*n2(3));
    % 
    % hs_tube = surf(Xt, Yt, Zt, ...
    %     'FaceAlpha', 0.10, 'EdgeColor', 'none', 'FaceColor', [0 0.4 1]);
    % 
    % % (optional) end-caps to make the tube extent obvious
    % % capTh = linspace(0, 2*pi, 64);
    % % C0 = ctr + t0_new*dstar_new + r_enclose_new*(cos(capTh)'*n1 + sin(capTh)'*n2);
    % % C1 = ctr + t1_new*dstar_new + r_enclose_new*(cos(capTh)'*n1 + sin(capTh)'*n2);
    % % patch('XData',C0(:,1),'YData',C0(:,2),'ZData',C0(:,3), ...
    % %       'FaceColor',[0 0.4 1],'FaceAlpha',0.10,'EdgeColor','none');
    % % patch('XData',C1(:,1),'YData',C1(:,2),'ZData',C1(:,3), ...
    % %       'FaceColor',[0 0.4 1],'FaceAlpha',0.10,'EdgeColor','none');
    % 
    % axis equal vis3d; grid on; box on; camlight headlight; lighting gouraud;
    % xlabel('X_{RAS} (mm)'); ylabel('Y_{RAS} (mm)'); zlabel('Z_{RAS} (mm)');
    % title('Planned trajectory (blue), fitted tube (blue translucent), grown blob (orange)');
    % view(135,20);

end