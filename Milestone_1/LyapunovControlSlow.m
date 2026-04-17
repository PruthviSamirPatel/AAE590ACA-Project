function [u, dbg] = LyapunovControlSlow(~, X, P)

r = X(1:3);
v = X(4:6);
r_norm = norm(r);

%% Current orbital elements
[a, ecc, inc, argp, raan, ~] = cart2kep(r(1), r(2), r(3), ...
                                         v(1), v(2), v(3), P.mu_nd);

%% Slow-variable errors
e_a       = (a - P.Target.a_nd) / P.Target.a_nd;
e_e       = ecc - P.Target.ecc;
e_i_deg   = wrapTo180_local(inc  - P.Target.inc);
e_raan_deg= wrapTo180_local(raan - P.Target.raan);
e_argp_deg= wrapTo180_local(argp - P.Target.argp);

e_i    = deg2rad(e_i_deg);
e_raan = deg2rad(e_raan_deg);
e_argp = deg2rad(e_argp_deg);

%% RTN frame
C_eci2rtn = eci2rtn_dcm(r, v);
C_rtn2eci = C_eci2rtn.';

%% Simple slow-variable control in RTN
% Radial mostly shapes e and omega
u_r = -P.k_e * e_e - P.k_argp * e_argp;

% Tangential mostly shapes a
u_t = -P.k_a * e_a;

% Normal mostly shapes i and RAAN
u_n = -P.k_i * e_i - P.k_raan * e_raan;

u_rtn = [u_r; u_t; u_n];

%% Simple penalty: outward radial push near Earth
if r_norm < P.r_warn
    sigma = (P.r_warn - r_norm) / (P.r_warn - P.r_min);
    sigma = max(0, sigma);
    u_pen_r = P.k_pen * sigma^2;
else
    u_pen_r = 0;
end

u_pen_rtn = [u_pen_r; 0; 0];

%% Total control in ECI
u_raw = C_rtn2eci * (u_rtn + u_pen_rtn);

%% Saturation
u_norm = norm(u_raw);
if u_norm > P.umax
    u = P.umax * u_raw / u_norm;
else
    u = u_raw;
end

%% Lyapunov-like diagnostic
Vslow = 0.5*(e_a^2 + e_e^2 + e_i^2 + e_raan^2 + e_argp^2);
Vpen  = 0.5*u_pen_r^2;
V     = Vslow + Vpen;

dbg.V = V;
dbg.penalty = Vpen;
dbg.e_a = e_a;
dbg.e_e = e_e;
dbg.e_i_deg = e_i_deg;
dbg.e_raan_deg = e_raan_deg;
dbg.e_argp_deg = e_argp_deg;
end