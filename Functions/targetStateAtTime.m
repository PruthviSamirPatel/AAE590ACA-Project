function [rT, vT, aT, M_deg] = targetStateAtTime(t, Target, mu_nd)

% Mean anomaly at time t
M_deg = Target.M + rad2deg(Target.n_nd * t);
M_deg = mod(M_deg, 360);

% Eccentric anomaly
E = kepler(M_deg, Target.ecc);

% True anomaly
ta = 2 * atan2( sqrt(1 + Target.ecc) * sind(E/2), ...
                sqrt(1 - Target.ecc) * cosd(E/2) );
ta = rad2deg(ta);

% Cartesian state in nondim units
[x, y, z, vx, vy, vz] = kep2cart(Target.a_nd, Target.ecc, Target.inc, ...
    Target.argp, Target.raan, ta, mu_nd);

rT = [x; y; z];
vT = [vx; vy; vz];
aT = -mu_nd * rT / norm(rT)^3;
end