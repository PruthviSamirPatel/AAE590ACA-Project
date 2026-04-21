clc; clear; close all

%% ===================== USER SETTINGS =====================
tFinal_days = 10;              % max time allowed for leg 2 transfer [days]

% Max control
uMax_km_s2 = 1e-3;             % km/s^2  (= 1 m/s^2)
                                % if you want 10 cm/s^2, this is correct
                                % if you want 1 cm/s^2, use 1e-5

% Earth avoidance / safety barrier settings
hSafe_km = 150;                % minimum allowed altitude above Earth surface
hOn_km   = 300;                % altitude where barrier starts turning on
kBarrier = 5e-5;               % barrier gain
epsBar   = 1e-6;               % avoids singularity in barrier

% Lyapunov gains
K = eye(3);
P = eye(3);

% Rendezvous event tolerances
rTol_km   = 1e-1;               % terminate when relative position < rTol_km
vTol_km_s = 1e-3;              % terminate when relative speed < vTol_km_s

%% ===================== LOAD DATA =====================
load("Earth_params.mat")
load("Orbital_Slot.mat")
load("Parking_Orbit.mat")
load("Intermediate_Orbit.mat")   % should contain InterAchieved from leg 1

rEarth = Earth.radius;          % km
mu     = Earth.mu;              % km^3/s^2
J2     = Earth.J2;

%% ===================== NONDIMENSIONALIZATION =====================
lStar = rEarth;
tStar = sqrt(rEarth^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;

mu_nd = 1;
Re_nd = rEarth / lStar;         % = 1 with this scaling
uMax  = uMax_km_s2 / aStar;     % nondimensionalized max acceleration

rSafe = (rEarth + hSafe_km) / lStar;
rOn   = (rEarth + hOn_km)   / lStar;

rTol_nd = rTol_km / lStar;
vTol_nd = vTol_km_s / vStar;

%% ===================== COAST TO RAAN MATCH =====================
[tMatch, raanI_match, raanT_match] = find_raan_match_time(Earth, InterAchieved, Target);

fprintf('RAAN match time = %.6f days\n', tMatch/(24*3600));
fprintf('Matched RAANs   = %.6f deg (inter) , %.6f deg (target)\n', raanI_match, raanT_match);

InterCoasted  = propagate_orbit_J2(Earth, InterAchieved, tMatch);
TargetCoasted = propagate_orbit_J2(Earth, Target,        tMatch);

%% ===================== INITIAL / TARGET ORBITS FOR LEG 2 =====================
% Controlled spacecraft starts from the coasted intermediary orbit
a0    = InterCoasted.a;
e0    = InterCoasted.ecc;
inc0  = InterCoasted.inc;     % deg
raan0 = InterCoasted.raan;    % deg
argp0 = InterCoasted.argp;    % deg
M0    = InterCoasted.M;       % deg

% Target spacecraft state after same coast time
aT    = TargetCoasted.a;
eT    = TargetCoasted.ecc;
incT  = TargetCoasted.inc;    % deg
raanT = TargetCoasted.raan;   % deg
argpT = TargetCoasted.argp;   % deg
MT    = TargetCoasted.M;      % deg

%% ===================== CONVERT MEAN ANOMALY TO TRUE ANOMALY =====================
% This assumes your kep2cart expects TRUE anomaly in degrees.
E0_deg = kepler(M0, e0);
ta0_deg = 2*atan2d(sqrt(1+e0)*sind(E0_deg/2), sqrt(1-e0)*cosd(E0_deg/2));
ta0_deg = mod(ta0_deg, 360);

ET_deg = kepler(MT, eT);
taT_deg = 2*atan2d(sqrt(1+eT)*sind(ET_deg/2), sqrt(1-eT)*cosd(ET_deg/2));
taT_deg = mod(taT_deg, 360);

%% ===================== INITIAL AND TARGET CARTESIAN STATES =====================
[x, y, z, vx, vy, vz] = kep2cart(a0, e0, inc0, argp0, raan0, ta0_deg, mu);
rC0_nd = [x; y; z] / lStar;
vC0_nd = [vx; vy; vz] / vStar;

[x, y, z, vx, vy, vz] = kep2cart(aT, eT, incT, argpT, raanT, taT_deg, mu);
rT0_nd = [x; y; z] / lStar;
vT0_nd = [vx; vy; vz] / vStar;

% Full state = [rC; vC; rT; vT]
X0 = [rC0_nd; vC0_nd; rT0_nd; vT0_nd];

%% ===================== TIME SPAN / EVENTS =====================
tSpan = [0, tFinal_days*24*3600/tStar];

opts = odeset('RelTol',1e-12, 'AbsTol',1e-12, ...
              'Events', @(t,X) combined_events_cart(t, X, rSafe, rTol_nd, vTol_nd));

%% ===================== INTEGRATE LEG 2 =====================
[t, X, te, Xe, ie] = ode45(@(t,X) cart_lyap_dyn(t, X, K, P, mu_nd, uMax, ...
                                                rSafe, rOn, kBarrier, epsBar, ...
                                                J2, Re_nd), ...
                           tSpan, X0, opts);

%% ===================== RECOVER HISTORIES =====================
rC_nd      = X(:,1:3);
vC_nd      = X(:,4:6);
rT_nd_hist = X(:,7:9);
vT_nd_hist = X(:,10:12);

rC_km = rC_nd * lStar;
vC_km_s = vC_nd * vStar;

rT_km = rT_nd_hist * lStar;
vT_km_s = vT_nd_hist * vStar;

dr_km   = (rC_nd - rT_nd_hist) * lStar;
dv_km_s = (vC_nd - vT_nd_hist) * vStar;

t_days = t * tStar / (24*3600);
t_sec  = t * tStar;

num_steps = length(t);
u_hist    = zeros(num_steps, 3);
u_norm    = zeros(num_steps, 1);
rmag_hist = zeros(num_steps, 1);
alt_hist  = zeros(num_steps, 1);

for ii = 1:num_steps
    r_i = X(ii,1:3)';

    rmag_hist(ii) = norm(r_i) * lStar;
    alt_hist(ii)  = rmag_hist(ii) - rEarth;

    [~, u_out] = cart_lyap_dyn_plot_helper(t(ii), X(ii,:)', K, P, mu_nd, ...
                                           uMax, rSafe, rOn, kBarrier, epsBar, ...
                                           J2, Re_nd);
    u_hist(ii,:) = u_out' * aStar;     % km/s^2
    u_norm(ii)   = norm(u_hist(ii,:));
end

%% ===================== DELTA-V =====================
DeltaV_km_s = trapz(t_sec, u_norm);
DeltaV_m_s  = 1000 * DeltaV_km_s;

fprintf('\nLeg 2 total Delta-V = %.6f m/s\n', DeltaV_m_s);
fprintf('Final relative position = %.6f km\n', norm(dr_km(end,:)));
fprintf('Final relative speed    = %.9f km/s\n', norm(dv_km_s(end,:)));

if ~isempty(ie)
    fprintf('Integration terminated by event ID = %d\n', ie(end));
    if ie(end) == 1
        fprintf('Reason: Earth keep-out boundary crossed.\n');
    elseif ie(end) == 2
        fprintf('Reason: Rendezvous tolerance achieved.\n');
    end
else
    fprintf('Integration ended at final time with no terminal event.\n');
end

%% ===================== SAVE FINAL STATE IF DESIRED =====================
Leg2Final.rC_km   = rC_km(end,:)';
Leg2Final.vC_km_s = vC_km_s(end,:)';
Leg2Final.rT_km   = rT_km(end,:)';
Leg2Final.vT_km_s = vT_km_s(end,:)';
Leg2Final.dr_km   = dr_km(end,:)';
Leg2Final.dv_km_s = dv_km_s(end,:)';
Leg2Final.DeltaV_m_s = DeltaV_m_s;
Leg2Final.t_days = t_days(end);

%% ===================== PLOTS =====================

% 3D trajectory plot
figure('Color','w');
hold on; grid on; axis equal; view(3);

[sx, sy, sz] = sphere(60);
topo = load('topo.mat');
surf(lStar*sx, lStar*sy, lStar*sz, 'FaceColor', 'texturemap', ...
     'CData', topo.topo, 'EdgeColor', 'none', 'HandleVisibility','off');

surf(rSafe*lStar*sx, rSafe*lStar*sy, rSafe*lStar*sz, ...
    'FaceAlpha', 0.08, 'EdgeColor', 'none', 'FaceColor', [1 0 0], ...
    'HandleVisibility','off');

plot3(rC_km(:,1), rC_km(:,2), rC_km(:,3), 'r', 'LineWidth', 1.5, ...
      'DisplayName', 'Controlled Spacecraft');
plot3(rT_km(:,1), rT_km(:,2), rT_km(:,3), 'b--', 'LineWidth', 1.5, ...
      'DisplayName', 'Target Spacecraft');

plot3(rC_km(1,1), rC_km(1,2), rC_km(1,3), 'ro', 'MarkerFaceColor','r', ...
      'DisplayName', 'Chaser Start');
plot3(rT_km(1,1), rT_km(1,2), rT_km(1,3), 'bo', 'MarkerFaceColor','b', ...
      'DisplayName', 'Target Start');

plot3(rC_km(end,1), rC_km(end,2), rC_km(end,3), 'rs', 'MarkerFaceColor','r', ...
      'DisplayName', 'Chaser End');
plot3(rT_km(end,1), rT_km(end,2), rT_km(end,3), 'bs', 'MarkerFaceColor','b', ...
      'DisplayName', 'Target End');

xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Leg 2 Cartesian Lyapunov Transfer with J2');
legend('show', 'Location', 'best');

% Relative position history
figure('Color','w');
plot(t_days, dr_km, 'LineWidth', 1.4);
grid on;
xlabel('t [days]');
ylabel('\Delta r [km]');
title('Relative Position History');
legend('\Deltax','\Deltay','\Deltaz','Location','best');

% Relative speed history
figure('Color','w');
plot(t_days, dv_km_s, 'LineWidth', 1.4); hold on;
plot(t_days, vecnorm(dv_km_s,2,2), 'k', 'LineWidth', 1.6);
grid on;
xlabel('t [days]');
ylabel('\Delta v [km/s]');
title('Relative Velocity History');
legend('\Deltav_x','\Deltav_y','\Deltav_z','||\Deltav||','Location','best');

% Control history
figure('Color','w');
plot(t_days, 1000*u_hist, 'LineWidth', 1.2); hold on;
plot(t_days, 1000*u_norm, 'k', 'LineWidth', 1.6);
grid on;
xlabel('t [days]');
ylabel('u [m/s^2]');
title('Control History');
legend('u_x','u_y','u_z','||u||','Location','best');

% Altitude history
figure('Color','w');
plot(t_days, alt_hist, 'LineWidth', 1.5); hold on;
yline(hSafe_km, 'r--', 'LineWidth', 1.5, 'DisplayName','Safe Altitude');
yline(hOn_km,   'k--', 'LineWidth', 1.2, 'DisplayName','Barrier On Altitude');
grid on;
xlabel('t [days]');
ylabel('Altitude [km]');
title('Controlled Spacecraft Altitude');
legend('Trajectory Altitude','Safe Altitude','Barrier On','Location','best');

% Norms of relative state
figure('Color','w');
plot(t_days, vecnorm(dr_km,2,2), 'LineWidth', 1.5); hold on;
plot(t_days, vecnorm(dv_km_s,2,2), 'LineWidth', 1.5);
grid on;
xlabel('t [days]');
ylabel('Norm');
title('Relative State Norms');
legend('||\Deltar|| [km]', '||\Deltav|| [km/s]', 'Location','best');

%% ========================= FUNCTIONS =========================

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

function [Xdot, u_dimless] = cart_lyap_dyn_plot_helper(~, X, K, P, mu, ...
                                                       uMax, rSafe, rOn, ...
                                                       kBarrier, epsBar, ...
                                                       J2, Re)

    rC = X(1:3);
    vC = X(4:6);
    rT = X(7:9);
    vT = X(10:12);

    dr = rC - rT;
    dv = vC - vT;

    g2B_C = -mu * rC / norm(rC)^3;
    g2B_T = -mu * rT / norm(rT)^3;

    aJ2_C = accelJ2_cart(rC, mu, J2, Re);
    aJ2_T = accelJ2_cart(rT, mu, J2, Re);

    gC = g2B_C + aJ2_C;
    gT = g2B_T + aJ2_T;
    delta_g = gC - gT;

    uLyap = -K*dr - P*dv - delta_g;

    rmag = norm(rC);
    d    = rmag - rSafe;
    dOn  = rOn - rSafe;

    uBarrier = [0;0;0];
    if d <= 0
        uBarrier = (kBarrier / epsBar) * (rC / norm(rC));
    elseif d < dOn
        uBarrier = kBarrier * (1/(d + epsBar) - 1/dOn) * (rC / norm(rC));
    end

    u = uLyap + uBarrier;

    if norm(u) > uMax
        u = uMax * u / norm(u);
    end

    u_dimless = u;

    rCdot = vC;
    vCdot = gC + u;

    rTdot = vT;
    vTdot = gT;

    Xdot = [rCdot; vCdot; rTdot; vTdot];
end

function aJ2 = accelJ2_cart(r, mu, J2, Re)
    x = r(1);
    y = r(2);
    z = r(3);

    r2 = x^2 + y^2 + z^2;
    r1 = sqrt(r2);
    z2 = z^2;

    factor = -(3/2) * J2 * mu * Re^2 / r1^5;
    common = 1 - 5*z2/r2;

    ax = factor * x * common;
    ay = factor * y * common;
    az = factor * z * (3 - 5*z2/r2);

    aJ2 = [ax; ay; az];
end

function [value, isterminal, direction] = combined_events_cart(~, X, rSafe, rTol, vTol)
    % Event 1: Earth keep-out
    rC = X(1:3);
    rmag = norm(rC);

    value1 = rmag - rSafe;

    % Event 2: Rendezvous achieved
    dr = X(1:3) - X(7:9);
    dv = X(4:6) - X(10:12);

    value2 = max(norm(dr) - rTol, norm(dv) - vTol);

    value = [value1; value2];
    isterminal = [1; 1];
    direction  = [-1; -1];
end

