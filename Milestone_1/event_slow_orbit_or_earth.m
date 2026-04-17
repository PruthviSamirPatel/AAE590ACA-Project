function [value, isterminal, direction] = event_slow_orbit_or_earth(~, X, P)

r = X(1:3);
v = X(4:6);

[a, ecc, inc, argp, raan, ~] = cart2kep(r(1), r(2), r(3), ...
                                         v(1), v(2), v(3), P.mu_nd);

e_a       = abs((a - P.Target.a_nd) / P.Target.a_nd);
e_e       = abs(ecc - P.Target.ecc);
e_i       = abs(deg2rad(wrapTo180_local(inc  - P.Target.inc)));
e_raan    = abs(deg2rad(wrapTo180_local(raan - P.Target.raan)));
e_argp    = abs(deg2rad(wrapTo180_local(argp - P.Target.argp)));

% Event 1: slow-variable convergence
tol_a    = 1e-3;
tol_e    = 1e-3;
tol_ang  = deg2rad(0.5);

value1 = max([e_a - tol_a, ...
              e_e - tol_e, ...
              e_i - tol_ang, ...
              e_raan - tol_ang, ...
              e_argp - tol_ang]);

% Event 2: Earth/path boundary
value2 = norm(r) - P.r_min;

value = [value1; value2];
isterminal = [1; 1];
direction = [-1; -1];
end