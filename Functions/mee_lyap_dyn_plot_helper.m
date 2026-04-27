function [Xdot, u_dimless] = mee_lyap_dyn_plot_helper(~, X, xslowT_nd, K, P, ...
                                                      mu, uMax, rSafe, rOn, ...
                                                      kBarrier, epsBar, J2, Re_nd)

    p = X(1);
    f = X(2);
    g = X(3);
    h = X(4);
    k = X(5);
    L = X(6);

    if p <= 0 || ~isfinite(p)
        Xdot = zeros(6,1);
        u_dimless = zeros(3,1);
        return
    end

    w  = 1 + f*cos(L) + g*sin(L);
    s2 = 1 + h^2 + k^2;
    r  = p / w;

    if w <= 1e-10 || r <= 0 || ~isfinite(r)
        Xdot = zeros(6,1);
        u_dimless = zeros(3,1);
        return
    end

    sq = sqrt(p/mu);

    f0 = [0; 0; 0; 0; 0; sqrt(mu/p^3)*w^2];

    B = zeros(6,3);

    B(1,:) = [0,            2*p/w*sq,                          0];

    B(2,:) = [sq*sin(L),    sq*((w+1)*cos(L) + f)/w, ...
                         -sq*(h*sin(L) - k*cos(L))/w];

    B(3,:) = [-sq*cos(L),   sq*((w+1)*sin(L) + g)/w, ...
                          sq*(h*cos(L) + k*sin(L))/w];

    B(4,:) = [0,            0,  sq*(s2/(2*w))*cos(L)];

    B(5,:) = [0,            0,  sq*(s2/(2*w))*sin(L)];

    B(6,:) = [0,            0,  sq*(h*sin(L) - k*cos(L))/w];

    Bslow = B(1:5,:);

    xslow = X(1:5);
    dx_slow = xslow - xslowT_nd;

    H = Bslow' * K' * K * Bslow + 1e-12*eye(3);
    gvec = Bslow' * K' * P * dx_slow;
    uLyap = -(H \ gvec);

    d   = r - rSafe;
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        uBarrier(1) = kBarrier / epsBar;
    elseif d < dOn
        uBarrier(1) = kBarrier * (1/(d + epsBar) - 1/dOn);
    end

    u = uLyap + uBarrier;

    unorm = norm(u);
    if unorm > uMax
        u = uMax * u / unorm;
    end

    u_dimless = u;
    [r_eci, v_eci] = mee2rv(X, mu);
    aJ2_eci = accelJ2(r_eci, mu, J2, Re_nd);
    aJ2_rtn = eci2rtn_accel(r_eci, v_eci, aJ2_eci);
    
    Xdot = f0 + B*(u + aJ2_rtn);
end