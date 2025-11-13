function C = point_at_s(model, E, D, C_bez, Cpoly, s_query)
% Return Nx3 points at arc positions s_query (mm from proximal)
switch lower(model)
    case 'line'
        L = max(norm(D-E),eps);
        t = min(max(s_query./L,0),1);
        C = E + t.*(D - E);
    case 'bezier'
        % sample densely, then interp by arc length
        [~, Cdense] = bezier_samples(E, C_bez, D, 1000);
        [s,~,L] = arc_length_param(Cdense);
        t = min(max(s_query./max(L,eps),0),1);
        C = [interp1(s,Cdense(:,1),t*L), ...
             interp1(s,Cdense(:,2),t*L), ...
             interp1(s,Cdense(:,3),t*L)];
    otherwise % 'polyline'
        [s,~,L] = arc_length_param(Cpoly);
        ss = min(max(s_query,0),L);
        C = [interp1(s,Cpoly(:,1),ss), ...
             interp1(s,Cpoly(:,2),ss), ...
             interp1(s,Cpoly(:,3),ss)];
end
end
function [s, seglen, Ltot] = arc_length_param(Cpts)
% Cpts: Mx3 points along centerline (ordered)
seglen = vecnorm(diff(Cpts,1,1),2,2);      % (M-1)x1
s = [0; cumsum(seglen)];                   % Mx1, mm
Ltot = s(end);
end