function [Xdot, u_dimless] = cart_lyap_dyn_plot_helper(~, X, K, P, mu, ...
                                                       uMax, rSafe, rOn, ...
                                                       kBarrier, epsBar, ...
                                                       J2, Re)

    rC = X(1:3);
    vC = X(4:6);
    rT = X(7:9);
    vT = X(10:12);

    dr = rC - rT;
    dv = vC - vT;

    g2B_C = -mu * rC / norm(rC)^3;
    g2B_T = -mu * rT / norm(rT)^3;

    aJ2_C = accelJ2_cart(rC, mu, J2, Re);
    aJ2_T = accelJ2_cart(rT, mu, J2, Re);

    gC = g2B_C + aJ2_C;
    gT = g2B_T + aJ2_T;
    delta_g = gC - gT;

    uLyap = -K*dr - P*dv - delta_g;

    rmag = norm(rC);
    d    = rmag - rSafe;
    dOn  = rOn - rSafe;

    uBarrier = [0;0;0];
    if d <= 0
        uBarrier = (kBarrier / epsBar) * (rC / norm(rC));
    elseif d < dOn
        uBarrier = kBarrier * (1/(d + epsBar) - 1/dOn) * (rC / norm(rC));
    end

    u = uLyap + uBarrier;

    if norm(u) > uMax
        u = uMax * u / norm(u);
    end

    u_dimless = u;

    rCdot = vC;
    vCdot = gC + u;

    rTdot = vT;
    vTdot = gT;

    Xdot = [rCdot; vCdot; rTdot; vTdot];
end