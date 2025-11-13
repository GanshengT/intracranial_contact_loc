function u = tan_at_s(ss, s, Cpts)
    % central difference in arc-length domain
    h = 0.5;
    s0 = max(s(1), ss - h);
    s1 = min(s(end), ss + h);
    p0 = [interp1(s,Cpts(:,1),s0), interp1(s,Cpts(:,2),s0), interp1(s,Cpts(:,3),s0)];
    p1 = [interp1(s,Cpts(:,1),s1), interp1(s,Cpts(:,2),s1), interp1(s,Cpts(:,3),s1)];
    u  = p1 - p0; n = norm(u); if n>0, u = u/n; else, u = [0 0 1]; end
end
% h
% Chooses a small step (in mm of arc-length) for a finite-difference.
% 	•	s0 = max(s(1), ss - h);  s1 = min(s(end), ss + h);
% Picks two nearby arc-lengths around ss, clamped to the valid range [s(1),
% s(end)] to avoid running off the curve near the ends.
% 	•	p0 = [...] interp1(s, Cpts(:,k), s0), p1 = [...] interp1(..., s1)
% Interpolates the 3D point on the centerline at those two arc-lengths
% (component-wise). This gives two nearby positions on the curve:
% \mathbf{p}_0=\mathbf{C}(s_0), \mathbf{p}_1=\mathbf{C}(s_1).
% 	•	u  = p1 - p0;  n = norm(u);  if n>0, u = u/n; else, u = [0 0 1];
% 	end
% Forms the central difference \mathbf{u} \propto
% \mathbf{C}(s_1)-\mathbf{C}(s_0) and normalizes it to unit length. If the
% difference degenerates (e.g., both points identical), it falls back to [0
% 0 1].