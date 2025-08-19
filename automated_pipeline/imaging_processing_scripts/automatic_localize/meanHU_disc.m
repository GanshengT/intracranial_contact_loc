function mhu = meanHU_disc(vol, Avox2ras, C, u, rad_mm, half_ax_mm)
    % Average HU over a small cylinder centered at RAS C,
    % axis along unit tangent u (RAS), radius rad_mm, half-axial half_ax_mm.
    C = C';
    u = u(:)/norm(u);
    % Orthonormal basis {u,v,w}
    [~,idx] = max(abs(u)); e = zeros(3,1); e(idx) = 1;
    v = e - (u'*e)*u; v = v/norm(v);
    w = cross(u,v);

    Nr=5; Nt=24; Na=3;
    rr = linspace(0, rad_mm, Nr);
    tt = linspace(0, 2*pi, Nt+1); tt(end)=[];
    aa = linspace(-half_ax_mm, +half_ax_mm, Na);

    P = zeros(Nr*Nt*Na,3); c=1;
    for a = aa
        for it = 1:Nt
            ct = cos(tt(it)); st = sin(tt(it));
            for ir = 1:Nr
                P(c,:) = C + a*u + rr(ir)*ct*v + rr(ir)*st*w;
                c = c+1;
            end
        end
    end
    % Sample nearest (robust to indexing)
    vals = sample_ct_nn(vol, Avox2ras, P);
    vals = vals(~isnan(vals));
    if isempty(vals), mhu = NaN; else, mhu = mean(vals); end
end