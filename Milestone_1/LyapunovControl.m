function [u, dbg] = LyapunovControl(t, X, P)

r = X(1:3);
v = X(4:6);

%% Target inertial state at time t
[rT, vT, aT, MT_deg] = targetStateAtTime(t, P.Target, P.mu_nd);

%% Current gravity
r_norm = norm(r);
aC = -P.mu_nd * r / r_norm^3;

%% Cartesian tracking errors
e_r = r - rT;
e_v = v - vT;

%% Base Lyapunov tracking term
u_track = -P.Kr * e_r - P.Kv * e_v + (aT - aC);

%% Current orbital elements
[a, ecc, inc, argp, raan, ta] = cart2kep(r(1), r(2), r(3), v(1), v(2), v(3), P.mu_nd);

if ~isreal(ecc) || ~isfinite(ecc) || ecc >= 1 || ecc < 0
    fprintf('\nInvalid eccentricity encountered in LyapunovControl:\n');
    fprintf('  t        = %.12f\n', t);
    fprintf('  ecc      = %.12e\n', ecc);
    fprintf('  a        = %.12e\n', a);
    fprintf('  ta       = %.12e\n', ta);
    fprintf('  r_norm   = %.12e\n', norm(r));
    fprintf('  v_norm   = %.12e\n', norm(v));
    disp('  r = '), disp(r.')
    disp('  v = '), disp(v.')
    error('LyapunovControl: non-elliptic or invalid state reached.');
end

% Mean anomaly from current state
M_deg = true2mean(ta, ecc);

% Angular errors
e_raan_deg = wrapTo180_local(raan - P.Target.raan);
e_M_deg    = wrapTo180_local(M_deg - MT_deg);

% Semimajor-axis relative error
e_a = (a - P.Target.a_nd) / P.Target.a_nd;

%% Orbital-element shaping in RTN frame
C_eci2rtn = eci2rtn_dcm(r, v);
C_rtn2eci = C_eci2rtn.';

u_rtn = zeros(3,1);

% Tangential thrust: mainly handles a and phasing
u_rtn(2) = -P.k_a * e_a - P.k_M * deg2rad(e_M_deg);

% Normal thrust: mainly handles RAAN correction
u_rtn(3) = -P.k_raan * deg2rad(e_raan_deg);

u_elem = C_rtn2eci * u_rtn;

%% Exponential path-constraint penalties
rhat = r / r_norm;

% Altitude barrier
% active mainly when r < r_warn, but remains smooth
delta_alt = P.r_warn - r_norm;
if delta_alt > 0
    u_barrier_alt = P.k_barrier * (exp(P.alpha_alt * delta_alt) - 1) * rhat;
    pen_alt = exp(P.alpha_alt * delta_alt) - 1;
else
    u_barrier_alt = zeros(3,1);
    pen_alt = 0;
end

% Perigee barrier
rp = a * (1 - ecc);     % nondimensional perigee radius
delta_rp = P.r_min - rp;
if delta_rp > 0
    u_barrier_rp = P.k_perigee * (exp(P.alpha_rp * delta_rp) - 1) * rhat;
    pen_rp = exp(P.alpha_rp * delta_rp) - 1;
else
    u_barrier_rp = zeros(3,1);
    pen_rp = 0;
end

u_barrier = u_barrier_alt + u_barrier_rp;

%% Total unsaturated control
u_raw = u_track + u_elem + u_barrier;

%% Saturation
u_norm = norm(u_raw);
if u_norm > P.umax
    u = P.umax * u_raw / u_norm;
else
    u = u_raw;
end

%% Lyapunov-like diagnostic
V_track = 0.5 * (e_r.' * P.Kr * e_r) + 0.5 * (e_v.' * e_v);
V_elem  = 0.5 * (e_a^2 + deg2rad(e_raan_deg)^2 + deg2rad(e_M_deg)^2);
V_pen   = 0.5 * (pen_alt^2 + pen_rp^2);
V = V_track + V_elem + V_pen;

dbg.V = V;
dbg.penalty = V_pen;

dbg.V = V;
dbg.penalty = V_pen;
dbg.e_a = e_a;
dbg.e_raan_deg = e_raan_deg;
dbg.e_M_deg = e_M_deg;
dbg.rT = rT;
dbg.vT = vT;
dbg.u_track = u_track;
dbg.u_elem = u_elem;
dbg.u_barrier = u_barrier;

end