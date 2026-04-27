clc; clear; close all

%% ===================== USER SETTINGS =====================
tFinal_days = 10;

% Max control
uMax_km_s2 = 1e-3;             % km/s^2

% Earth avoidance
hSafe_km = 150;
hOn_km   = 300;
kBarrier = 5e-5;
epsBar   = 1e-6;

% Full MEE Lyapunov gains
K = eye(6);
P = eye(6);

% Make the fast angle target less aggressive if needed
K(6,6) = 0.25;
P(6,6) = 0.25;

% Terminal tolerances
meeTol = [1e-5; 1e-5; 1e-5; 1e-5; 1e-5; deg2rad(0.1)];

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

%% ===================== COAST TO RAAN MATCH =====================
[tMatch, raanI_match, raanT_match] = find_raan_match_time(Earth, InterAchieved, Target);

fprintf('RAAN match time = %.6f days\n', tMatch/(24*3600));
fprintf('Matched RAANs   = %.6f deg inter, %.6f deg target\n', raanI_match, raanT_match);

InterCoasted  = propagate_orbit_J2(Earth, InterAchieved, tMatch);
TargetCoasted = propagate_orbit_J2(Earth, Target,        tMatch);

%% ===================== INITIAL / TARGET MEE STATES =====================
% Chaser starts from coasted intermediary orbit
aC    = InterCoasted.a;
eC    = InterCoasted.ecc;
incC  = deg2rad(InterCoasted.inc);
raanC = deg2rad(InterCoasted.raan);
argpC = deg2rad(InterCoasted.argp);
MC    = InterCoasted.M;

% Target starts from coasted target slot
aT    = TargetCoasted.a;
eT    = TargetCoasted.ecc;
incT  = deg2rad(TargetCoasted.inc);
raanT = deg2rad(TargetCoasted.raan);
argpT = deg2rad(TargetCoasted.argp);
MT    = TargetCoasted.M;

% Mean anomaly to true anomaly
EC_deg = kepler(MC, eC);
taC = deg2rad(mod(2*atan2d(sqrt(1+eC)*sind(EC_deg/2), ...
                           sqrt(1-eC)*cosd(EC_deg/2)), 360));

ET_deg = kepler(MT, eT);
taT = deg2rad(mod(2*atan2d(sqrt(1+eT)*sind(ET_deg/2), ...
                           sqrt(1-eT)*cosd(ET_deg/2)), 360));

% Dimensional MEE
meeC_dim = kep2mee_local(aC, eC, incC, raanC, argpC, taC);
meeT_dim = kep2mee_local(aT, eT, incT, raanT, argpT, taT);

% Nondimensional MEE: only p is length-scaled
meeC0 = [meeC_dim(1)/lStar; meeC_dim(2:6)];
meeT0 = [meeT_dim(1)/lStar; meeT_dim(2:6)];

% State is [controlled MEE; target MEE]
X0 = [meeC0; meeT0];

%% ===================== TIME SPAN / EVENTS =====================
tSpan = [0, tFinal_days*24*3600/tStar];

opts = odeset('RelTol',1e-12, 'AbsTol',1e-12, ...
    'Events', @(t,X) combined_events_full_mee(t, X, rSafe, meeTol));

%% ===================== INTEGRATE LEG 2 =====================
[t, X, te, Xe, ie] = ode45(@(t,X) full_mee_leg2_dyn(t, X, K, P, mu_nd, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd), tSpan, X0, opts);

%% ===================== RECOVER HISTORIES =====================
num_steps = length(t);
t_days = t*tStar/(24*3600);
t_sec  = t*tStar;

meeC_hist = X(:,1:6);
meeT_hist = X(:,7:12);

rC_km = zeros(num_steps,3);
vC_km_s = zeros(num_steps,3);
rT_km = zeros(num_steps,3);
vT_km_s = zeros(num_steps,3);

u_hist = zeros(num_steps,3);
u_norm = zeros(num_steps,1);
alt_hist = zeros(num_steps,1);

a_hist = zeros(num_steps,1);
e_hist = zeros(num_steps,1);
inc_hist = zeros(num_steps,1);
raan_hist = zeros(num_steps,1);
argp_hist = zeros(num_steps,1);
ta_hist = zeros(num_steps,1);

mee_err = zeros(num_steps,6);

for i = 1:num_steps
    meeC_i = meeC_hist(i,:)';
    meeT_i = meeT_hist(i,:)';

    [rC_nd, vC_nd] = mee2rv_local(meeC_i, mu_nd);
    [rT_nd, vT_nd] = mee2rv_local(meeT_i, mu_nd);

    rC_km(i,:) = rC_nd'*lStar;
    vC_km_s(i,:) = vC_nd'*vStar;
    rT_km(i,:) = rT_nd'*lStar;
    vT_km_s(i,:) = vT_nd'*vStar;

    alt_hist(i) = norm(rC_nd)*lStar - rEarth;

    meeC_dim_i = [meeC_i(1)*lStar; meeC_i(2:6)];
    [a_i,e_i,inc_i,raan_i,argp_i,ta_i] = mee2kep_local(meeC_dim_i);

    a_hist(i) = a_i;
    e_hist(i) = e_i;
    inc_hist(i) = inc_i;
    raan_hist(i) = raan_i;
    argp_hist(i) = argp_i;
    ta_hist(i) = ta_i;

    mee_err(i,:) = full_mee_error(meeC_i, meeT_i)';

    [~, u_out] = full_mee_leg2_plot_helper(t(i), X(i,:)', K, P, mu_nd, uMax, ...
        rSafe, rOn, kBarrier, epsBar, J2, Re_nd);

    u_hist(i,:) = u_out'*aStar;
    u_norm(i) = norm(u_hist(i,:));
end

dr_km = rC_km - rT_km;
dv_km_s = vC_km_s - vT_km_s;

DeltaV_km_s = trapz(t_sec, u_norm);
DeltaV_m_s = 1000*DeltaV_km_s;

fprintf('\nLeg 2 total Delta-V = %.6f m/s\n', DeltaV_m_s);
fprintf('Final relative position = %.6f km\n', norm(dr_km(end,:)));
fprintf('Final relative speed    = %.9f km/s\n', norm(dv_km_s(end,:)));
fprintf('Final full MEE error norm = %.6e\n', norm(mee_err(end,:)));

if ~isempty(ie)
    fprintf('Integration terminated by event ID = %d\n', ie(end));
    if ie(end) == 1
        fprintf('Reason: Earth keep-out boundary crossed.\n');
    elseif ie(end) == 2
        fprintf('Reason: Full MEE tolerance achieved.\n');
    end
else
    fprintf('Integration ended at final time with no terminal event.\n');
end

%% ===================== SAVE FINAL STATE =====================
Leg2Final.rC_km = rC_km(end,:)';
Leg2Final.vC_km_s = vC_km_s(end,:)';
Leg2Final.rT_km = rT_km(end,:)';
Leg2Final.vT_km_s = vT_km_s(end,:)';
Leg2Final.dr_km = dr_km(end,:)';
Leg2Final.dv_km_s = dv_km_s(end,:)';
Leg2Final.meeC_nd = meeC_hist(end,:)';
Leg2Final.meeT_nd = meeT_hist(end,:)';
Leg2Final.mee_err = mee_err(end,:)';
Leg2Final.DeltaV_m_s = DeltaV_m_s;
Leg2Final.t_days = t_days(end);

%% ===================== PLOTS =====================
figure('Color','w');
hold on; grid on; axis equal; view(3);

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
    'DisplayName','Target Spacecraft');

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
title('Leg 2 Full-MEE Lyapunov Transfer with J2');
legend('show','Location','best');

figure('Color','w');
plot(t_days, mee_err, 'LineWidth',1.4);
grid on;
xlabel('t [days]');
ylabel('MEE error');
title('Full MEE Error History');
legend('\Deltap','\Deltaf','\Deltag','\Deltah','\Deltak','\DeltaL','Location','best');

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
title(sprintf('Control History — \\DeltaV = %.3f km/s', DeltaV_km_s));
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

figure('Color','w');
subplot(3,2,1)
plot(t_days, a_hist, 'LineWidth',1.5); grid on
xlabel('t [days]'); ylabel('a [km]')

subplot(3,2,2)
plot(t_days, e_hist, 'LineWidth',1.5); grid on
xlabel('t [days]'); ylabel('e')

subplot(3,2,3)
plot(t_days, rad2deg(inc_hist), 'LineWidth',1.5); grid on
xlabel('t [days]'); ylabel('i [deg]')

subplot(3,2,4)
plot(t_days, rad2deg(raan_hist), 'LineWidth',1.5); grid on
xlabel('t [days]'); ylabel('\Omega [deg]')

subplot(3,2,5)
plot(t_days, rad2deg(argp_hist), 'LineWidth',1.5); grid on
xlabel('t [days]'); ylabel('\omega [deg]')

subplot(3,2,6)
plot(t_days, rad2deg(ta_hist), 'LineWidth',1.5); grid on
xlabel('t [days]'); ylabel('\nu [deg]')

sgtitle('Recovered Classical Elements of Controlled Spacecraft')

%% ========================= FUNCTIONS =========================

function Xdot = full_mee_leg2_dyn(~, X, K, P, mu, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd)

    meeC = X(1:6);
    meeT = X(7:12);

    [fC, BC, aJ2C_rtn, rC] = mee_natural_terms(meeC, mu, J2, Re_nd);
    [fT, BT, aJ2T_rtn, ~]  = mee_natural_terms(meeT, mu, J2, Re_nd);

    xerr = full_mee_error(meeC, meeT);

    natRel = (fC + BC*aJ2C_rtn) - (fT + BT*aJ2T_rtn);

    H = BC' * K' * K * BC + 1e-12*eye(3);
    rhs = BC' * K' * K * (natRel + P*xerr);

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

function [Xdot, u_dimless] = full_mee_leg2_plot_helper(~, X, K, P, mu, uMax, ...
    rSafe, rOn, kBarrier, epsBar, J2, Re_nd)

    meeC = X(1:6);
    meeT = X(7:12);

    [fC, BC, aJ2C_rtn, rC] = mee_natural_terms(meeC, mu, J2, Re_nd);
    [fT, BT, aJ2T_rtn, ~]  = mee_natural_terms(meeT, mu, J2, Re_nd);

    xerr = full_mee_error(meeC, meeT);

    natRel = (fC + BC*aJ2C_rtn) - (fT + BT*aJ2T_rtn);

    H = BC' * K' * K * BC + 1e-12*eye(3);
    rhs = BC' * K' * K * (natRel + P*xerr);

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

function err = full_mee_error(meeC, meeT)
    err = meeC - meeT;
    err(6) = wrapToPiLocal(meeC(6) - meeT(6));
end

function [value, isterminal, direction] = combined_events_full_mee(~, X, rSafe, meeTol)

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

    err = abs(full_mee_error(meeC, meeT));
    value2 = max(err./meeTol) - 1;

    value = [value1; value2];
    isterminal = [1; 1];
    direction = [-1; -1];
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

function ang = wrapToPiLocal(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end

function ang = wrapTo2PiLocal(ang)
    ang = mod(ang, 2*pi);
end