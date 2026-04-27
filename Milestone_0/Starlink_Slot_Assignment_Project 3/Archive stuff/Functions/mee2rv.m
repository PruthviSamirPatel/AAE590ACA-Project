function [r_eci, v_eci] = mee2rv(X, mu)
% Convert modified equinoctial elements to inertial Cartesian state
% X = [p; f; g; h; k; L]

    p = X(1);
    f = X(2);
    g = X(3);
    h = X(4);
    k = X(5);
    L = X(6);

    cL = cos(L);
    sL = sin(L);

    w = 1 + f*cL + g*sL;
    s2 = 1 + h^2 + k^2;
    r = p / w;

    % Position and velocity in equinoctial frame
    x_pf = r * cL;
    y_pf = r * sL;

    vx_pf = -sqrt(mu/p) * (g + sL);
    vy_pf =  sqrt(mu/p) * (f + cL);

    % Basis vectors
    fhat = 1/s2 * [1 - h^2 + k^2;
                   2*h*k;
                  -2*h];

    ghat = 1/s2 * [2*h*k;
                   1 + h^2 - k^2;
                   2*k];

    % Inertial position and velocity
    r_eci = x_pf*fhat + y_pf*ghat;
    v_eci = vx_pf*fhat + vy_pf*ghat;
end