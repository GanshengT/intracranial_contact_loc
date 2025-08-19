function [v1,v2] = orthobasis_row(u)
    u = reshape(u,1,3); u = u/max(norm(u),eps);
    [~,idx] = max(abs(u)); e = zeros(1,3); e(idx)=1;
    v1 = cross(u,e); if norm(v1)<eps, v1 = cross(u,[1 0 0]); end
    v1 = v1/max(norm(v1),eps);
    v2 = cross(u,v1); v2 = v2/max(norm(v2),eps);
end