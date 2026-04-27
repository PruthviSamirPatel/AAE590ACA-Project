function Xdot = cart_lyap_dyn(~, X, K, P, mu, uMax, ...
                              rSafe, rOn, kBarrier, epsBar, J2, Re)

    % Unpack state
    rC = X(1:3);
    vC = X(4:6);
    rT = X(7:9);
    vT = X(10:12);

    % Relative state
    dr = rC - rT;
    dv = vC - vT;

    % Two-body gravity
    g2B_C = -mu * rC / norm(rC)^3;
    g2B_T = -mu * rT / norm(rT)^3;

    % J2 accelerations
    aJ2_C = accelJ2_cart(rC, mu, J2, Re);
    aJ2_T = accelJ2_cart(rT, mu, J2, Re);

    % Total natural accelerations
    gC = g2B_C + aJ2_C;
    gT = g2B_T + aJ2_T;
    delta_g = gC - gT;

    % Lyapunov control
    uLyap = -K*dr - P*dv - delta_g;

    % Earth barrier
    rmag = norm(rC);
    d    = rmag - rSafe;
    dOn  = rOn - rSafe;

    uBarrier = [0;0;0];
    if d <= 0
        uBarrier = (kBarrier / epsBar) * (rC / norm(rC));
    elseif d < dOn
        uBarrier = kBarrier * (1/(d + epsBar) - 1/dOn) * (rC / norm(rC));
    end

    % Commanded control
    u = uLyap + uBarrier;

    % Saturation on control only
    if norm(u) > uMax
        u = uMax * u / norm(u);
    end

    % Dynamics
    rCdot = vC;
    vCdot = gC + u;

    rTdot = vT;
    vTdot = gT;

    Xdot = [rCdot; vCdot; rTdot; vTdot];
end