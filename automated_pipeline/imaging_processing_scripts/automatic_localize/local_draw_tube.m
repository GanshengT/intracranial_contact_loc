function local_draw_tube(ctr, d, t0, t1, r_mm)
    EPS = 1e-9;
    d = d(:).'/max(norm(d),EPS);
    % orthonormal frame {d, n1, n2}
    tmp = [1 0 0]; if abs(dot(tmp,d))>0.9, tmp = [0 1 0]; end
    n1 = tmp - dot(tmp,d)*d; n1 = n1 / max(norm(n1),EPS);
    n2 = cross(d, n1);       n2 = n2 / max(norm(n2),EPS);
    L  = max(t1 - t0, 1e-3);
    nTau = max(24, ceil(L/3)); nTh = 48;
    tau   = linspace(t0, t1, nTau);
    theta = linspace(0, 2*pi, nTh);
    [TH, TAU] = meshgrid(theta, tau);
    Xt = ctr(1) + TAU.*d(1) + r_mm*(cos(TH)*n1(1) + sin(TH)*n2(1));
    Yt = ctr(2) + TAU.*d(2) + r_mm*(cos(TH)*n1(2) + sin(TH)*n2(2));
    Zt = ctr(3) + TAU.*d(3) + r_mm*(cos(TH)*n1(3) + sin(TH)*n2(3));
    surf(Xt, Yt, Zt, 'FaceAlpha',0.10, 'EdgeColor','none'); 
    capTh = linspace(0,2*pi,64);
    C0 = ctr + t0*d + r_mm*(cos(capTh)'*n1 + sin(capTh)'*n2);
    C1 = ctr + t1*d + r_mm*(cos(capTh)'*n1 + sin(capTh)'*n2);
    patch('XData',C0(:,1),'YData',C0(:,2),'ZData',C0(:,3), 'FaceAlpha',0.10,'EdgeColor','none');
    patch('XData',C1(:,1),'YData',C1(:,2),'ZData',C1(:,3), 'FaceAlpha',0.10,'EdgeColor','none');
end