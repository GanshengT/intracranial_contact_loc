function [V,F] = cylinder_at_pose(C,u,L,r,nTheta,nSeg)
    [v1,v2,u] = axis_frame(u);
    % segment endpoints
    E = C - 0.5*L*u;
    D = C + 0.5*L*u;
    [V,F] = tube_mesh_segment(E, D, r, nTheta, nSeg); % from earlier helper
end

function [v1,v2,u] = axis_frame(u)
    u = u(:).'; u = u / max(norm(u),eps);
    ref = [1 0 0]; if abs(dot(ref,u))>0.9, ref=[0 1 0]; end
    v1 = cross(u, ref); v1 = v1 / max(norm(v1),eps);
    v2 = cross(u, v1);
end