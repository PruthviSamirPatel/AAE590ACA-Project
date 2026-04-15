clc; clear; close all;

fprintf('\n============================================================\n');
fprintf(' Lyapunov Reference Trajectory Generation with Path Penalty\n');
fprintf('============================================================\n');

%% Load environment
load("Earth_params.mat")
rEarth = Earth.radius;          % km
mu     = Earth.mu;              % km^3/s^2

%% Load target and parking orbit
load("Orbital_Slot.mat")
load("Parking_Orbit.mat")

%% Nondimensionalization
lStar = rEarth;
tStar = sqrt(rEarth^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;
mu_nd = 1;

fprintf('\nEnvironment loaded:\n');
fprintf('  Earth radius          = %.6f km\n', rEarth);
fprintf('  Earth mu              = %.6f km^3/s^2\n', mu);
fprintf('  lStar                 = %.6f km\n', lStar);
fprintf('  tStar                 = %.6f s\n', tStar);
fprintf('  vStar                 = %.6f km/s\n', vStar);
fprintf('  aStar                 = %.9f km/s^2\n', aStar);

%% Control authority
umax_dim = 1e-5;                 % km/s^2
umax     = umax_dim / aStar;     % nondimensional

fprintf('\nControl authority:\n');
fprintf('  umax_dim              = %.6e km/s^2\n', umax_dim);
fprintf('  umax_nd               = %.6f\n', umax);

%% Initial state from parking orbit
E0 = kepler(Parking.M, Parking.ecc);
ta0 = 2 * atan2( sqrt(1 + Parking.ecc) * sind(E0/2), ...
                 sqrt(1 - Parking.ecc) * cosd(E0/2) );
ta0 = rad2deg(ta0);

[x0, y0, z0, vx0, vy0, vz0] = kep2cart(Parking.a, Parking.ecc, Parking.inc, ...
    Parking.argp, Parking.raan, ta0, mu);

r0_dim = [x0; y0; z0];
v0_dim = [vx0; vy0; vz0];

r0 = r0_dim / lStar;
v0 = v0_dim / vStar;
X0 = [r0; v0];

fprintf('\nInitial parking orbit:\n');
fprintf('  a0                    = %.6f km\n', Parking.a);
fprintf('  e0                    = %.8f\n', Parking.ecc);
fprintf('  i0                    = %.6f deg\n', Parking.inc);
fprintf('  RAAN0                 = %.6f deg\n', Parking.raan);
fprintf('  omega0                = %.6f deg\n', Parking.argp);
fprintf('  M0                    = %.6f deg\n', Parking.M);
fprintf('  TA0                   = %.6f deg\n', ta0);

fprintf('\nTarget orbit at initial epoch:\n');
fprintf('  aT                    = %.6f km\n', Target.a);
fprintf('  eT                    = %.8f\n', Target.ecc);
fprintf('  iT                    = %.6f deg\n', Target.inc);
fprintf('  RAANT                 = %.6f deg\n', Target.raan);
fprintf('  omegaT                = %.6f deg\n', Target.argp);
fprintf('  MT(0)                 = %.6f deg\n', Target.M);

%% Guidance parameters
P.mu_nd   = mu_nd;
P.umax    = umax;
P.Target  = Target;

% Cartesian Lyapunov gains
P.Kr = diag([0.05, 0.05, 0.05]);
P.Kv = diag([0.20, 0.20, 0.20]);

% Orbital-element shaping gains (RTN acceleration commands)
P.k_a    = 0.01;
P.k_M    = 0.005;
P.k_raan = 0.005;

% Earth avoidance / path penalty parameters
h_min_km      = 150;                     % minimum acceptable altitude
h_warn_km     = 400;                     % warning altitude
P.r_min       = (rEarth + h_min_km) / lStar;
P.r_warn      = (rEarth + h_warn_km) / lStar;

% Exponential penalty gains
P.k_barrier   = 8.0;                     % altitude barrier strength
P.k_perigee   = 4.0;                     % perigee barrier strength
P.alpha_alt   = 60.0;                    % exponential steepness for altitude
P.alpha_rp    = 80.0;                    % exponential steepness for perigee

fprintf('\nGuidance parameters:\n');
fprintf('  Kr diag               = [%.3f %.3f %.3f]\n', diag(P.Kr));
fprintf('  Kv diag               = [%.3f %.3f %.3f]\n', diag(P.Kv));
fprintf('  k_a                   = %.6f\n', P.k_a);
fprintf('  k_M                   = %.6f\n', P.k_M);
fprintf('  k_raan                = %.6f\n', P.k_raan);
fprintf('  h_min                 = %.3f km\n', h_min_km);
fprintf('  h_warn                = %.3f km\n', h_warn_km);
fprintf('  r_min_nd              = %.9f\n', P.r_min);
fprintf('  r_warn_nd             = %.9f\n', P.r_warn);
fprintf('  k_barrier             = %.6f\n', P.k_barrier);
fprintf('  k_perigee             = %.6f\n', P.k_perigee);
fprintf('  alpha_alt             = %.6f\n', P.alpha_alt);
fprintf('  alpha_rp              = %.6f\n', P.alpha_rp);

%% Integration settings
t0 = 0;
tf = 30;     % nondimensional final horizon

opts = odeset('RelTol',1e-10, 'AbsTol',1e-10, ...
              'Events', @(t,X) event_target_or_earth(t,X,P));

fprintf('\nIntegration settings:\n');
fprintf('  t0                    = %.6f nondim\n', t0);
fprintf('  tf                    = %.6f nondim\n', tf);
fprintf('  tf                    = %.6f s\n', tf*tStar);

%% Propagate Lyapunov-guided trajectory
fprintf('\nPropagating Lyapunov-guided reference trajectory...\n');

[tSol, XSol, tEvent, XEvent, iEvent] = ode45(@(t,X) EoM_LyapunovPenalty(t, X, P), ...
    [t0 tf], X0, opts);

fprintf('Propagation complete.\n');

%% Recover control history and diagnostics
N = length(tSol);

uHist        = zeros(N,3);
VHist        = zeros(N,1);
penHist      = zeros(N,1);
aErrHist     = zeros(N,1);
raanErrHist  = zeros(N,1);
MErrHist     = zeros(N,1);
altHist_km   = zeros(N,1);
rpHist_km    = zeros(N,1);
aHist_km     = zeros(N,1);
eHist        = zeros(N,1);
iHist_deg    = zeros(N,1);
raanHist_deg = zeros(N,1);
argpHist_deg = zeros(N,1);
taHist_deg   = zeros(N,1);
MHist_deg    = zeros(N,1);

for k = 1:N
    Xk = XSol(k,:).';
    rk = Xk(1:3);
    vk = Xk(4:6);

    [u, dbg] = LyapunovControl(tSol(k), Xk, P);

    uHist(k,:)       = u.';
    VHist(k)         = dbg.V;
    penHist(k)       = dbg.penalty;
    aErrHist(k)      = dbg.e_a;
    raanErrHist(k)   = dbg.e_raan_deg;
    MErrHist(k)      = dbg.e_M_deg;

    [a_k, e_k, i_k, argp_k, raan_k, ta_k] = cart2kep(rk(1), rk(2), rk(3), ...
                                                      vk(1), vk(2), vk(3), mu_nd);

    aHist_km(k)     = a_k * lStar;
    eHist(k)        = e_k;
    iHist_deg(k)    = i_k;
    raanHist_deg(k) = raan_k;
    argpHist_deg(k) = argp_k;
    taHist_deg(k)   = ta_k;
    MHist_deg(k)    = true2mean(ta_k, e_k);

    altHist_km(k)   = norm(rk)*lStar - rEarth;
    rpHist_km(k)    = a_k*(1 - e_k)*lStar - rEarth;
end

%% Final state and final orbit
tf_actual = tSol(end);
Xf = XSol(end,:).';
rf = Xf(1:3);
vf = Xf(4:6);

[a_f_nd, e_f, i_f, argp_f, raan_f, ta_f] = cart2kep(rf(1), rf(2), rf(3), ...
                                                     vf(1), vf(2), vf(3), mu_nd);
a_f_km = a_f_nd * lStar;
M_f_deg = true2mean(ta_f, e_f);

[rT_f, vT_f, ~, MT_f_deg] = targetStateAtTime(tf_actual, Target, mu_nd);
[aT_f_nd, eT_f, iT_f, argpT_f, raanT_f, taT_f] = cart2kep(rT_f(1), rT_f(2), rT_f(3), ...
                                                          vT_f(1), vT_f(2), vT_f(3), mu_nd);
aT_f_km = aT_f_nd * lStar;

% Final tracking errors
pos_err_nd = norm(rf - rT_f);
vel_err_nd = norm(vf - vT_f);
pos_err_km = pos_err_nd * lStar;
vel_err_kms = vel_err_nd * vStar;

e_a_final = (a_f_nd - aT_f_nd)/aT_f_nd;
e_raan_final = wrapTo180_local(raan_f - raanT_f);
e_M_final    = wrapTo180_local(M_f_deg - MT_f_deg);

% Min altitude/perigee
minAlt_km = min(altHist_km);
minRpAlt_km = min(rpHist_km);

% Control stats
uMag = vecnorm(uHist,2,2);
uMaxUsed_nd = max(uMag);
uMaxUsed_dim = uMaxUsed_nd * aStar;
uMean_nd = mean(uMag);
uMean_dim = uMean_nd * aStar;

% Delta-v estimate
dv_nd = trapz(tSol, uMag);
dv_dim = dv_nd * vStar;

%% Event reporting
fprintf('\n============================================================\n');
fprintf(' EVENT REPORT\n');
fprintf('============================================================\n');

if isempty(iEvent)
    fprintf('No terminal event occurred. Integration reached tf.\n');
else
    fprintf('Integration terminated by event index: %d\n', iEvent(end));
    switch iEvent(end)
        case 1
            fprintf('Meaning: target tolerance reached.\n');
        case 2
            fprintf('Meaning: Earth/path boundary reached.\n');
        otherwise
            fprintf('Meaning: unrecognized event index.\n');
    end
    fprintf('Event time              = %.6f nondim\n', tEvent(end));
    fprintf('Event time              = %.6f s\n', tEvent(end)*tStar);
end

%% Print important results
fprintf('\n============================================================\n');
fprintf(' IMPORTANT TRANSFER RESULTS\n');
fprintf('============================================================\n');
fprintf('Final integration time   = %.6f nondim\n', tf_actual);
fprintf('Final integration time   = %.6f s\n', tf_actual*tStar);
fprintf('Final integration time   = %.6f hr\n', tf_actual*tStar/3600);
fprintf('Minimum altitude         = %.6f km\n', minAlt_km);
fprintf('Minimum perigee altitude = %.6f km\n', minRpAlt_km);
fprintf('Max control used         = %.6e nondim\n', uMaxUsed_nd);
fprintf('Max control used         = %.6e km/s^2\n', uMaxUsed_dim);
fprintf('Mean control magnitude   = %.6e nondim\n', uMean_nd);
fprintf('Mean control magnitude   = %.6e km/s^2\n', uMean_dim);
fprintf('Approx. delta-v          = %.6f km/s\n', dv_dim);
fprintf('Final position error     = %.6f km\n', pos_err_km);
fprintf('Final velocity error     = %.6f km/s\n', vel_err_kms);
fprintf('Final rel. a error       = %.6e\n', e_a_final);
fprintf('Final RAAN error         = %.6f deg\n', e_raan_final);
fprintf('Final mean anomaly error = %.6f deg\n', e_M_final);
fprintf('Final V                  = %.6e\n', VHist(end));
fprintf('Final penalty            = %.6e\n', penHist(end));

%% Print final reference orbit
fprintf('\n============================================================\n');
fprintf(' FINAL REFERENCE ORBIT\n');
fprintf('============================================================\n');
fprintf('Reference orbit at t_f:\n');
fprintf('  a_ref                 = %.6f km\n', a_f_km);
fprintf('  e_ref                 = %.8f\n', e_f);
fprintf('  i_ref                 = %.6f deg\n', i_f);
fprintf('  RAAN_ref              = %.6f deg\n', raan_f);
fprintf('  omega_ref             = %.6f deg\n', argp_f);
fprintf('  TA_ref                = %.6f deg\n', ta_f);
fprintf('  M_ref                 = %.6f deg\n', M_f_deg);

fprintf('\nTarget orbit at t_f:\n');
fprintf('  a_tgt(tf)             = %.6f km\n', aT_f_km);
fprintf('  e_tgt(tf)             = %.8f\n', eT_f);
fprintf('  i_tgt(tf)             = %.6f deg\n', iT_f);
fprintf('  RAAN_tgt(tf)          = %.6f deg\n', raanT_f);
fprintf('  omega_tgt(tf)         = %.6f deg\n', argpT_f);
fprintf('  TA_tgt(tf)            = %.6f deg\n', taT_f);
fprintf('  M_tgt(tf)             = %.6f deg\n', MT_f_deg);

%% Plot 3D trajectory
figure;
plot3(XSol(:,1), XSol(:,2), XSol(:,3), 'b', 'LineWidth', 1.5); hold on;
[sx, sy, sz] = sphere(60);
surf(sx, sy, sz, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot3(XSol(1,1), XSol(1,2), XSol(1,3), 'go', 'MarkerFaceColor','g');
plot3(XSol(end,1), XSol(end,2), XSol(end,3), 'ro', 'MarkerFaceColor','r');
axis equal; grid on;
xlabel('x / R_E'); ylabel('y / R_E'); zlabel('z / R_E');
title('Lyapunov Reference Trajectory with Exponential Path Penalty');
legend('Trajectory','Earth','Start','Final');

%% Plot control
figure;
plot(tSol, uHist(:,1), 'LineWidth', 1.2); hold on;
plot(tSol, uHist(:,2), 'LineWidth', 1.2);
plot(tSol, uHist(:,3), 'LineWidth', 1.2);
plot(tSol, uMag, '--k', 'LineWidth', 1.3);
grid on;
xlabel('t / t^*');
ylabel('u (nondim)');
legend('u_x','u_y','u_z','||u||');
title('Control History');

%% Plot Lyapunov and penalty
figure;
plot(tSol, VHist, 'LineWidth', 1.5); hold on;
plot(tSol, penHist, '--', 'LineWidth', 1.5);
grid on;
xlabel('t / t^*');
ylabel('Value');
legend('V','Penalty');
title('Lyapunov Function and Exponential Path Penalty');

%% Plot orbital errors
figure;
subplot(3,1,1);
plot(tSol, aErrHist, 'LineWidth', 1.2); grid on;
ylabel('\Delta a / a_T');
title('Orbital-Element Tracking Errors');

subplot(3,1,2);
plot(tSol, raanErrHist, 'LineWidth', 1.2); grid on;
ylabel('\Delta \Omega [deg]');

subplot(3,1,3);
plot(tSol, MErrHist, 'LineWidth', 1.2); grid on;
ylabel('\Delta M [deg]');
xlabel('t / t^*');

%% Plot safety metrics
figure;
subplot(2,1,1);
plot(tSol, altHist_km, 'LineWidth', 1.3); hold on;
yline(h_min_km, '--r', 'h_{min}');
yline(h_warn_km, '--m', 'h_{warn}');
grid on;
ylabel('Altitude [km]');
title('Altitude History');

subplot(2,1,2);
plot(tSol, rpHist_km, 'LineWidth', 1.3); hold on;
yline(h_min_km, '--r', 'h_{min}');
grid on;
ylabel('Perigee altitude [km]');
xlabel('t / t^*');
title('Perigee Altitude History');

%% Save reference trajectory for later linearization
RefTraj.t              = tSol;
RefTraj.X              = XSol;
RefTraj.u              = uHist;
RefTraj.uMag           = uMag;
RefTraj.V              = VHist;
RefTraj.penalty        = penHist;
RefTraj.alt_km         = altHist_km;
RefTraj.rpAlt_km       = rpHist_km;
RefTraj.a_km           = aHist_km;
RefTraj.ecc            = eHist;
RefTraj.inc_deg        = iHist_deg;
RefTraj.raan_deg       = raanHist_deg;
RefTraj.argp_deg       = argpHist_deg;
RefTraj.ta_deg         = taHist_deg;
RefTraj.M_deg          = MHist_deg;
RefTraj.final_state    = Xf;
RefTraj.final_time_nd  = tf_actual;
RefTraj.final_time_s   = tf_actual*tStar;
RefTraj.final_orbit.a_km    = a_f_km;
RefTraj.final_orbit.ecc     = e_f;
RefTraj.final_orbit.inc_deg = i_f;
RefTraj.final_orbit.raan_deg = raan_f;
RefTraj.final_orbit.argp_deg = argp_f;
RefTraj.final_orbit.ta_deg   = ta_f;
RefTraj.final_orbit.M_deg    = M_f_deg;
RefTraj.target_final.a_km    = aT_f_km;
RefTraj.target_final.ecc     = eT_f;
RefTraj.target_final.inc_deg = iT_f;
RefTraj.target_final.raan_deg = raanT_f;
RefTraj.target_final.argp_deg = argpT_f;
RefTraj.target_final.ta_deg   = taT_f;
RefTraj.target_final.M_deg    = MT_f_deg;
RefTraj.pos_err_km      = pos_err_km;
RefTraj.vel_err_kms     = vel_err_kms;
RefTraj.e_a_final       = e_a_final;
RefTraj.e_raan_final    = e_raan_final;
RefTraj.e_M_final       = e_M_final;
RefTraj.dv_kms          = dv_dim;
RefTraj.minAlt_km       = minAlt_km;
RefTraj.minRpAlt_km     = minRpAlt_km;
RefTraj.event_index     = [];
RefTraj.event_time_nd   = [];
RefTraj.event_time_s    = [];

if ~isempty(iEvent)
    RefTraj.event_index   = iEvent(end);
    RefTraj.event_time_nd = tEvent(end);
    RefTraj.event_time_s  = tEvent(end)*tStar;
end

save('Reference_Trajectory_Lyapunov.mat', 'RefTraj');

fprintf('\nReference trajectory saved to: Reference_Trajectory_Lyapunov.mat\n');
fprintf('Done.\n\n');