function OrbitOut = propagate_orbit_J2(Earth, Orbit0, t)
% Propagate classical orbital elements under secular J2 effects
%
% INPUTS:
%   Earth  - struct with fields:
%            .mu     [km^3/s^2]
%            .J2     [-]
%            .radius [km]
%
%   Orbit0 - struct with fields:
%            .a      [km]
%            .ecc    [-]
%            .inc    [deg]
%            .raan   [deg]
%            .argp   [deg]
%            .M      [deg]
%
%   t      - propagation time [s]
%
% OUTPUT:
%   OrbitOut - propagated orbit struct with same fields plus drift rates

    % Unpack
    mu  = Earth.mu;
    J2  = Earth.J2;
    Re  = Earth.radius;

    a    = Orbit0.a;
    ecc  = Orbit0.ecc;
    inc  = Orbit0.inc;
    raan = Orbit0.raan;
    argp = Orbit0.argp;
    M    = Orbit0.M;

    % Basic quantities
    p = a*(1 - ecc^2);
    n = sqrt(mu/a^3);   % rad/s
    inc_rad = deg2rad(inc);

    % Secular J2 rates
    raanDot = -1.5 * J2 * (Re/p)^2 * n * cos(inc_rad);   % rad/s

    argpDot = 0.75 * J2 * (Re/p)^2 * n * (5*cos(inc_rad)^2 - 1); % rad/s

    MDot = n + 0.75 * J2 * (Re/p)^2 * n * sqrt(1 - ecc^2) * ...
              (3*cos(inc_rad)^2 - 1); % rad/s

    % Propagate
    OrbitOut = Orbit0;
    OrbitOut.a    = a;
    OrbitOut.ecc  = ecc;
    OrbitOut.inc  = inc;
    OrbitOut.raan = mod(raan + rad2deg(raanDot*t), 360);
    OrbitOut.argp = mod(argp + rad2deg(argpDot*t), 360);
    OrbitOut.M    = mod(M    + rad2deg(MDot*t),   360);
    OrbitOut.n    = n;

    % Save rates too
    OrbitOut.raanDot = raanDot;   % rad/s
    OrbitOut.argpDot = argpDot;   % rad/s
    OrbitOut.MDot    = MDot;      % rad/s
end