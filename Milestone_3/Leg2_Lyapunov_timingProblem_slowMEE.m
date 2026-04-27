clc; clear; close all

%% ===================== USER SETTINGS =====================
tFinal_days = 10;              % max allowed controlled Leg 2 time [days]

% Extra coast search after RAAN match
extraCoast_orbits_search = 2;  % search over this many target orbits
N_extra_grid = 41;             % grid points used to find fzero bracket

% Max control
uMax_km_s2 = 1e-3;             % km/s^2

% Earth avoidance
hSafe_km = 150;
hOn_km   = 300;
kBarrier = 5e-5;
epsBar   = 1e-6;

% Slow-MEE gains
K = eye(5);
P = eye(5);

% Slow-MEE terminal tolerances
slowTol = [1e-5; 1e-6; 1e-6; 1e-6; 1e-6];

%% ===================== LOAD DATA =====================
load("Earth_params.mat")
load("Orbital_Slot.mat")
load("Parking_Orbit.mat")
load("Intermediate_Orbit.mat")   % contains InterAchieved

rEarth = Earth.radius;
mu     = Earth.mu;
J2     = Earth.J2;

%% ===================== NONDIMENSIONALIZATION =====================
lStar = rEarth;
tStar = sqrt(rEarth^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;

mu_nd = 1;
Re_nd = rEarth/lStar;
uMax  = uMax_km_s2/aStar;

rSafe = (rEarth + hSafe_km)/lStar;
rOn   = (rEarth + hOn_km)/lStar;

tFinal_nd = tFinal_days*24*3600/tStar;

%% ===================== BASE RAAN MATCH COAST =====================
[tMatch, raanI_match, raanT_match] = find_raan_match_time(Earth, InterAchieved, Target);

fprintf('Base RAAN match time = %.6f days\n', tMatch/(24*3600));
fprintf('Matched RAANs = %.6f deg inter, %.6f deg target\n', ...
    raanI_match, raanT_match);

TargetAtMatch = propagate_orbit_J2(Earth, Target, tMatch);
nTarget = sqrt(mu/TargetAtMatch.a^3);
Ttarget = 2*pi/nTarget;

extra_min = 0;
extra_max = extraCoast_orbits_search*Ttarget;

%% ===================== OUTER FSOLVE FOR IDEAL EXTRA COAST =====================
fprintf('\nSolving for ideal extra coast time using fsolve...\n')

% Initial guess: half of one target orbit after RAAN match
extraGuess = 0.9*Ttarget;

phaseFun_fsolve = @(dtExtra) phase_error_for_fsolve(dtExtra, ...
    tMatch, Ttarget, Earth, InterAchieved, Target, ...
    lStar, tStar, vStar, aStar, mu_nd, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd, ...
    K, P, slowTol, tFinal_nd);

fsolveOpts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'MaxIterations',30, ...
    'MaxFunctionEvaluations',100);

extraCoast_best = fsolve(phaseFun_fsolve, extraGuess, fsolveOpts);

% Wrap to a physically reasonable equivalent coast in [0, Ttarget)
extraCoast_best = mod(extraCoast_best, Ttarget);

%% ===================== FINAL RUN WITH IDEAL EXTRA COAST =====================
Result = run_leg2_slowMEE_for_extra_coast(extraCoast_best, ...
    tMatch, Earth, InterAchieved, Target, ...
    lStar, tStar, vStar, aStar, mu_nd, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd, ...
    K, P, slowTol, tFinal_nd);

%% ===================== FINAL RUN WITH IDEAL EXTRA COAST =====================
Result = run_leg2_slowMEE_for_extra_coast(extraCoast_best, ...
    tMatch, Earth, InterAchieved, Target, ...
    lStar, tStar, vStar, aStar, mu_nd, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd, ...
    K, P, slowTol, tFinal_nd);

fprintf('\n===================== FINAL RESULT =====================\n')
fprintf('Base RAAN match coast = %.6f days\n', tMatch/(24*3600));
fprintf('Extra coast           = %.6f days\n', extraCoast_best/(24*3600));
fprintf('Total pre-transfer coast = %.6f days\n', (tMatch + extraCoast_best)/(24*3600));
fprintf('Controlled transfer time = %.6f days\n', Result.tTransfer_days);
fprintf('Final phase error = %.9f deg\n', rad2deg(Result.phaseErr_rad));
fprintf('Final mean anomaly error = %.9f deg\n', rad2deg(Result.Merr_rad));
fprintf('Final relative position = %.6f km\n', norm(Result.dr_km(end,:)));
fprintf('Final relative speed    = %.9f km/s\n', norm(Result.dv_km_s(end,:)));
fprintf('Leg 2 Delta-V = %.6f m/s\n', Result.DeltaV_m_s);

%% ===================== SAVE RESULT =====================
Leg2Final = Result;
save('Leg2_SlowMEE_Timed_Result.mat', 'Leg2Final')

%% ===================== PLOTS =====================
t_days = Result.t_days;
rC_km = Result.rC_km;
rT_km = Result.rT_km;
dr_km = Result.dr_km;
dv_km_s = Result.dv_km_s;
u_hist = Result.u_hist;
u_norm = Result.u_norm;
alt_hist = Result.alt_hist;
mee_err = Result.mee_err;

figure('Color','w');
hold on; grid on; axis equal; view(3)

[sx, sy, sz] = sphere(60);
topo = load('topo.mat');

surf(lStar*sx, lStar*sy, lStar*sz, 'FaceColor','texturemap', ...
    'CData',topo.topo, 'EdgeColor','none', 'HandleVisibility','off');

surf(rSafe*lStar*sx, rSafe*lStar*sy, rSafe*lStar*sz, ...
    'FaceAlpha',0.08, 'EdgeColor','none', 'FaceColor',[1 0 0], ...
    'HandleVisibility','off');

plot3(rC_km(:,1), rC_km(:,2), rC_km(:,3), 'r', 'LineWidth',1.5, ...
    'DisplayName','Controlled Spacecraft');

plot3(rT_km(:,1), rT_km(:,2), rT_km(:,3), 'b--', 'LineWidth',1.5, ...
    'DisplayName','Target Slot');

plot3(rC_km(1,1), rC_km(1,2), rC_km(1,3), 'ro', 'MarkerFaceColor','r', ...
    'DisplayName','Chaser Start');

plot3(rT_km(1,1), rT_km(1,2), rT_km(1,3), 'bo', 'MarkerFaceColor','b', ...
    'DisplayName','Target Start');

plot3(rC_km(end,1), rC_km(end,2), rC_km(end,3), 'rs', 'MarkerFaceColor','r', ...
    'DisplayName','Chaser End');

plot3(rT_km(end,1), rT_km(end,2), rT_km(end,3), 'bs', 'MarkerFaceColor','b', ...
    'DisplayName','Target End');

xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Leg 2 Slow-MEE Transfer with Timed Extra Coast');
legend('show','Location','best');

figure('Color','w');
plot(t_days, mee_err, 'LineWidth',1.4);
grid on;
xlabel('t [days]');
ylabel('Slow MEE Error');
title('Slow-MEE Error History');
legend('\Deltap','\Deltaf','\Deltag','\Deltah','\Deltak','Location','best');

figure('Color','w');
plot(t_days, dr_km, 'LineWidth',1.4);
grid on;
xlabel('t [days]');
ylabel('\Delta r [km]');
title('Relative Position History');
legend('\Deltax','\Deltay','\Deltaz','Location','best');

figure('Color','w');
plot(t_days, dv_km_s, 'LineWidth',1.4); hold on;
plot(t_days, vecnorm(dv_km_s,2,2), 'k', 'LineWidth',1.6);
grid on;
xlabel('t [days]');
ylabel('\Delta v [km/s]');
title('Relative Velocity History');
legend('\Deltav_x','\Deltav_y','\Deltav_z','||\Deltav||','Location','best');

figure('Color','w');
plot(t_days, 1000*u_hist, 'LineWidth',1.2); hold on;
plot(t_days, 1000*u_norm, 'k', 'LineWidth',1.6);
grid on;
xlabel('t [days]');
ylabel('u [m/s^2]');
title(sprintf('Control History — \\DeltaV = %.3f km/s', Result.DeltaV_m_s/1000));
legend('u_R','u_T','u_N','||u||','Location','best');

figure('Color','w');
plot(t_days, alt_hist, 'LineWidth',1.5); hold on;
yline(hSafe_km, 'r--', 'LineWidth',1.5);
yline(hOn_km, 'k--', 'LineWidth',1.2);
grid on;
xlabel('t [days]');
ylabel('Altitude [km]');
title('Controlled Spacecraft Altitude');
legend('Trajectory Altitude','Safe Altitude','Barrier On','Location','best');

%% ========================================================================
%% ============================= FUNCTIONS ================================
%% ========================================================================

function Result = run_leg2_slowMEE_for_extra_coast(dtExtra, ...
    tMatch, Earth, InterAchieved, Target, ...
    lStar, tStar, vStar, aStar, mu_nd, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd, ...
    K, P, slowTol, tFinal_nd)

    mu = Earth.mu;
    rEarth = Earth.radius;

    tPre = tMatch + dtExtra;

    InterStart  = propagate_orbit_J2(Earth, InterAchieved, tPre);
    TargetStart = propagate_orbit_J2(Earth, Target,        tPre);

    meeC0 = oe_struct_to_mee_nd(InterStart, mu, lStar);
    meeT0 = oe_struct_to_mee_nd(TargetStart, mu, lStar);

    X0 = [meeC0; meeT0];

    opts = odeset('RelTol',1e-12, 'AbsTol',1e-12, ...
        'Events', @(t,X) combined_events_slowMEE(t, X, rSafe, slowTol));

    [t, X, te, Xe, ie] = ode45(@(t,X) slowMEE_leg2_dyn(t, X, K, P, mu_nd, uMax, ...
        rSafe, rOn, kBarrier, epsBar, J2, Re_nd), [0, tFinal_nd], X0, opts);

    num_steps = length(t);

    meeC_hist = X(:,1:6);
    meeT_hist = X(:,7:12);

    rC_km = zeros(num_steps,3);
    vC_km_s = zeros(num_steps,3);
    rT_km = zeros(num_steps,3);
    vT_km_s = zeros(num_steps,3);
    u_hist = zeros(num_steps,3);
    u_norm = zeros(num_steps,1);
    alt_hist = zeros(num_steps,1);
    mee_err = zeros(num_steps,5);

    for i = 1:num_steps
        meeC = meeC_hist(i,:)';
        meeT = meeT_hist(i,:)';

        [rC_nd, vC_nd] = mee2rv_local(meeC, mu_nd);
        [rT_nd, vT_nd] = mee2rv_local(meeT, mu_nd);

        rC_km(i,:) = rC_nd'*lStar;
        vC_km_s(i,:) = vC_nd'*vStar;
        rT_km(i,:) = rT_nd'*lStar;
        vT_km_s(i,:) = vT_nd'*vStar;

        alt_hist(i) = norm(rC_nd)*lStar - rEarth;

        mee_err(i,:) = slow_mee_error(meeC, meeT)';

        [~, u_out] = slowMEE_leg2_plot_helper(t(i), X(i,:)', K, P, mu_nd, uMax, ...
            rSafe, rOn, kBarrier, epsBar, J2, Re_nd);

        u_hist(i,:) = u_out'*aStar;
        u_norm(i) = norm(u_hist(i,:));
    end

    dr_km = rC_km - rT_km;
    dv_km_s = vC_km_s - vT_km_s;

    t_sec = t*tStar;
    t_days = t_sec/(24*3600);

    DeltaV_km_s = trapz(t_sec, u_norm);
    DeltaV_m_s = 1000*DeltaV_km_s;

    meeC_final = meeC_hist(end,:)';
    meeT_final = meeT_hist(end,:)';

    phaseErr_rad = wrapToPiLocal(meeC_final(6) - meeT_final(6));

    meeC_dim = [meeC_final(1)*lStar; meeC_final(2:6)];
    meeT_dim = [meeT_final(1)*lStar; meeT_final(2:6)];

    [~, eC, ~, ~, ~, taC] = mee2kep_local(meeC_dim);
    [~, eT, ~, ~, ~, taT] = mee2kep_local(meeT_dim);

    MC = true_to_mean_anomaly(taC, eC);
    MT = true_to_mean_anomaly(taT, eT);

    Merr_rad = wrapToPiLocal(MC - MT);

    Result.dtExtra = dtExtra;
    Result.tPre = tPre;
    Result.t = t;
    Result.t_days = t_days;
    Result.tTransfer_days = t_days(end);

    Result.X = X;
    Result.meeC_hist = meeC_hist;
    Result.meeT_hist = meeT_hist;
    Result.mee_err = mee_err;

    Result.rC_km = rC_km;
    Result.vC_km_s = vC_km_s;
    Result.rT_km = rT_km;
    Result.vT_km_s = vT_km_s;
    Result.dr_km = dr_km;
    Result.dv_km_s = dv_km_s;

    Result.u_hist = u_hist;
    Result.u_norm = u_norm;
    Result.alt_hist = alt_hist;

    Result.phaseErr_rad = phaseErr_rad;
    Result.Merr_rad = Merr_rad;
    Result.DeltaV_m_s = DeltaV_m_s;

    Result.te = te;
    Result.Xe = Xe;
    Result.ie = ie;
end

function Xdot = slowMEE_leg2_dyn(~, X, K, P, mu, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd)

    meeC = X(1:6);
    meeT = X(7:12);

    [fC, BC, aJ2C_rtn, rC] = mee_natural_terms(meeC, mu, J2, Re_nd);
    [fT, BT, aJ2T_rtn, ~]  = mee_natural_terms(meeT, mu, J2, Re_nd);

    Bslow = BC(1:5,:);

    xerr = slow_mee_error(meeC, meeT);

    natRelSlow = (fC(1:5) + BC(1:5,:)*aJ2C_rtn) - ...
                 (fT(1:5) + BT(1:5,:)*aJ2T_rtn);

    H = Bslow' * K' * K * Bslow + 1e-12*eye(3);
    rhs = Bslow' * K' * K * (natRelSlow + P*xerr);

    uLyap = -(H\rhs);

    rmag = norm(rC);
    d = rmag - rSafe;
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        uBarrier(1) = kBarrier/epsBar;
    elseif d < dOn
        uBarrier(1) = kBarrier*(1/(d + epsBar) - 1/dOn);
    end

    u = uLyap + uBarrier;

    if norm(u) > uMax
        u = uMax*u/norm(u);
    end

    meeCdot = fC + BC*(u + aJ2C_rtn);
    meeTdot = fT + BT*aJ2T_rtn;

    Xdot = [meeCdot; meeTdot];
end

function [Xdot, u_dimless] = slowMEE_leg2_plot_helper(~, X, K, P, mu, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd)

    meeC = X(1:6);
    meeT = X(7:12);

    [fC, BC, aJ2C_rtn, rC] = mee_natural_terms(meeC, mu, J2, Re_nd);
    [fT, BT, aJ2T_rtn, ~]  = mee_natural_terms(meeT, mu, J2, Re_nd);

    Bslow = BC(1:5,:);

    xerr = slow_mee_error(meeC, meeT);

    natRelSlow = (fC(1:5) + BC(1:5,:)*aJ2C_rtn) - ...
                 (fT(1:5) + BT(1:5,:)*aJ2T_rtn);

    H = Bslow' * K' * K * Bslow + 1e-12*eye(3);
    rhs = Bslow' * K' * K * (natRelSlow + P*xerr);

    uLyap = -(H\rhs);

    rmag = norm(rC);
    d = rmag - rSafe;
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        uBarrier(1) = kBarrier/epsBar;
    elseif d < dOn
        uBarrier(1) = kBarrier*(1/(d + epsBar) - 1/dOn);
    end

    u = uLyap + uBarrier;

    if norm(u) > uMax
        u = uMax*u/norm(u);
    end

    u_dimless = u;

    meeCdot = fC + BC*(u + aJ2C_rtn);
    meeTdot = fT + BT*aJ2T_rtn;

    Xdot = [meeCdot; meeTdot];
end

function err = slow_mee_error(meeC, meeT)
    err = meeC(1:5) - meeT(1:5);
end

function [value, isterminal, direction] = combined_events_slowMEE(~, X, rSafe, slowTol)

    meeC = X(1:6);
    meeT = X(7:12);

    p = meeC(1);
    f = meeC(2);
    g = meeC(3);
    L = meeC(6);

    w = 1 + f*cos(L) + g*sin(L);

    if p <= 0 || w <= 0
        value1 = -1;
    else
        r = p/w;
        value1 = r - rSafe;
    end

    err = abs(slow_mee_error(meeC, meeT));
    value2 = max(err./slowTol) - 1;

    value = [value1; value2];
    isterminal = [1; 1];
    direction = [-1; -1];
end

function mee_nd = oe_struct_to_mee_nd(Orb, mu, lStar)

    a = Orb.a;
    e = Orb.ecc;
    inc = deg2rad(Orb.inc);
    raan = deg2rad(Orb.raan);
    argp = deg2rad(Orb.argp);
    M_deg = Orb.M;

    E_deg = kepler(M_deg, e);

    ta_deg = mod(2*atan2d(sqrt(1+e)*sind(E_deg/2), ...
                          sqrt(1-e)*cosd(E_deg/2)), 360);

    ta = deg2rad(ta_deg);

    mee_dim = kep2mee_local(a, e, inc, raan, argp, ta);
    mee_nd = [mee_dim(1)/lStar; mee_dim(2:6)];
end

function [f0, B, aJ2_rtn, r_eci] = mee_natural_terms(mee, mu, J2, Re_nd)

    p = mee(1);
    f = mee(2);
    g = mee(3);
    h = mee(4);
    k = mee(5);
    L = mee(6);

    w = 1 + f*cos(L) + g*sin(L);
    s2 = 1 + h^2 + k^2;

    if p <= 0 || w <= 1e-12
        f0 = zeros(6,1);
        B = zeros(6,3);
        aJ2_rtn = zeros(3,1);
        r_eci = zeros(3,1);
        return
    end

    sq = sqrt(p/mu);

    f0 = [0; 0; 0; 0; 0; sqrt(mu/p^3)*w^2];

    B = zeros(6,3);

    B(1,:) = [0, 2*p/w*sq, 0];

    B(2,:) = [sq*sin(L), ...
              sq*((w+1)*cos(L) + f)/w, ...
             -sq*(h*sin(L) - k*cos(L))/w];

    B(3,:) = [-sq*cos(L), ...
               sq*((w+1)*sin(L) + g)/w, ...
               sq*(h*cos(L) + k*sin(L))/w];

    B(4,:) = [0, 0, sq*(s2/(2*w))*cos(L)];

    B(5,:) = [0, 0, sq*(s2/(2*w))*sin(L)];

    B(6,:) = [0, 0, sq*(h*sin(L) - k*cos(L))/w];

    [r_eci, v_eci] = mee2rv_local(mee, mu);
    aJ2_eci = accelJ2_cart_local(r_eci, mu, J2, Re_nd);
    aJ2_rtn = eci2rtn_accel_local(r_eci, v_eci, aJ2_eci);
end

function mee = kep2mee_local(a, e, inc, raan, argp, ta)

    p = a*(1 - e^2);
    f = e*cos(argp + raan);
    g = e*sin(argp + raan);
    h = tan(inc/2)*cos(raan);
    k = tan(inc/2)*sin(raan);
    L = wrapTo2PiLocal(raan + argp + ta);

    mee = [p; f; g; h; k; L];
end

function [a,e,inc,raan,argp,ta] = mee2kep_local(mee)

    p = mee(1);
    f = mee(2);
    g = mee(3);
    h = mee(4);
    k = mee(5);
    L = mee(6);

    e = sqrt(f^2 + g^2);
    a = p/(1 - e^2);

    raan = atan2(k,h);
    inc = 2*atan(sqrt(h^2 + k^2));

    lonPer = atan2(g,f);
    argp = wrapTo2PiLocal(lonPer - raan);
    ta = wrapTo2PiLocal(L - lonPer);

    raan = wrapTo2PiLocal(raan);
end

function [r_eci, v_eci] = mee2rv_local(mee, mu)

    [a,e,inc,raan,argp,ta] = mee2kep_local(mee);

    p = a*(1 - e^2);

    r_pf = p/(1 + e*cos(ta))*[cos(ta); sin(ta); 0];
    v_pf = sqrt(mu/p)*[-sin(ta); e + cos(ta); 0];

    R3_W = [ cos(raan), -sin(raan), 0;
             sin(raan),  cos(raan), 0;
             0,          0,         1];

    R1_i = [1, 0,        0;
            0, cos(inc), -sin(inc);
            0, sin(inc),  cos(inc)];

    R3_w = [ cos(argp), -sin(argp), 0;
             sin(argp),  cos(argp), 0;
             0,          0,         1];

    Q = R3_W*R1_i*R3_w;

    r_eci = Q*r_pf;
    v_eci = Q*v_pf;
end

function aJ2 = accelJ2_cart_local(r, mu, J2, Re)

    x = r(1);
    y = r(2);
    z = r(3);

    r2 = x^2 + y^2 + z^2;
    r1 = sqrt(r2);
    z2 = z^2;

    factor = -(3/2)*J2*mu*Re^2/r1^5;
    common = 1 - 5*z2/r2;

    ax = factor*x*common;
    ay = factor*y*common;
    az = factor*z*(3 - 5*z2/r2);

    aJ2 = [ax; ay; az];
end

function a_rtn = eci2rtn_accel_local(r, v, a_eci)

    rhat = r/norm(r);
    hvec = cross(r,v);
    nhat = hvec/norm(hvec);
    that = cross(nhat,rhat);

    C = [rhat'; that'; nhat'];

    a_rtn = C*a_eci;
end

function M = true_to_mean_anomaly(ta, e)

    E = 2*atan2(sqrt(1-e)*sin(ta/2), ...
                sqrt(1+e)*cos(ta/2));

    E = wrapTo2PiLocal(E);
    M = E - e*sin(E);
    M = wrapTo2PiLocal(M);
end

function ang = wrapToPiLocal(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end

function ang = wrapTo2PiLocal(ang)
    ang = mod(ang, 2*pi);
end

function phaseErr = phase_error_for_fsolve(dtExtra, ...
    tMatch, Ttarget, Earth, InterAchieved, Target, ...
    lStar, tStar, vStar, aStar, mu_nd, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd, ...
    K, P, slowTol, tFinal_nd)

    % Keep fsolve from testing negative coast times
    dtExtra = mod(dtExtra, Ttarget);

    Result = run_leg2_slowMEE_for_extra_coast(dtExtra, ...
        tMatch, Earth, InterAchieved, Target, ...
        lStar, tStar, vStar, aStar, mu_nd, uMax, ...
        rSafe, rOn, kBarrier, epsBar, J2, Re_nd, ...
        K, P, slowTol, tFinal_nd);

    phaseErr = Result.phaseErr_rad;

    fprintf('fsolve trial: extra coast = %.9f days, phase error = %.9f deg\n', ...
        dtExtra/(24*3600), rad2deg(phaseErr));
end