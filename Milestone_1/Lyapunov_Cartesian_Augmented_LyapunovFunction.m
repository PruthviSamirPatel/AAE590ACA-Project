clc; clear; close all

%% Load environment
load("Earth_params.mat")
rEarth = Earth.radius;          % km
mu     = Earth.mu;              % km^3/s^2

%% Load target and parking orbit
load("Orbital_Slot.mat")
load("Parking_Orbit.mat")

%% Nondimensionalization
lStar = rEarth;
tStar = sqrt(lStar^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;
mu_nd = 1;

%% Max control
uMax_km_s2 = 1e-3;              % km/s^2
uMax = uMax_km_s2/aStar;        % nondimensional

%% Earth avoidance / augmented exponential penalty settings
hSafe_km = 100;                 % hard minimum altitude
hOn_km   = 200;                 % penalty becomes noticeable around here

rSafe = (rEarth + hSafe_km)/lStar;
rOn   = (rEarth + hOn_km)/lStar;

% Constraint:
% g = rMin - ||rC|| <= 0
rMin = rSafe;

% Makes penalty active before reaching rSafe
epsPenalty = rOn - rSafe;

% Exponential penalty: P(g) = exp(k*(g + eps))
wPenalty = 1.0;
kPenalty = 50.0;

%% Initial orbit
a0    = Parking.a;
e0    = Parking.ecc;
inc0  = Parking.inc;
raan0 = Parking.raan;
argp0 = Parking.argp;
M0    = Parking.M;

%% Target orbit
aT    = Target.a;
eT    = Target.ecc;
incT  = Target.inc;
raanT = Target.raan;
argpT = Target.argp;
MT    = Target.M;

%% Convert mean anomaly to true anomaly
E0_deg = kepler(M0, e0);
ta0_deg = 2*atan2d(sqrt(1+e0)*sind(E0_deg/2), ...
                   sqrt(1-e0)*cosd(E0_deg/2));
ta0_deg = mod(ta0_deg, 360);

ET_deg = kepler(MT, eT);
taT_deg = 2*atan2d(sqrt(1+eT)*sind(ET_deg/2), ...
                   sqrt(1-eT)*cosd(ET_deg/2));
taT_deg = mod(taT_deg, 360);

%% Initial and target Cartesian states
[x, y, z, vx, vy, vz] = kep2cart(a0, e0, inc0, argp0, raan0, ta0_deg, mu);
rC0_nd = [x; y; z]/lStar;
vC0_nd = [vx; vy; vz]/vStar;

[x, y, z, vx, vy, vz] = kep2cart(aT, eT, incT, argpT, raanT, taT_deg, mu);
rT0_nd = [x; y; z]/lStar;
vT0_nd = [vx; vy; vz]/vStar;

% Full state: [rC; vC; rT; vT]
X0 = [rC0_nd; vC0_nd; rT0_nd; vT0_nd];

%% Gains
K = eye(3);
P = eye(3);

%% Time span
tFinal_days = 10;
tSpan = [0, tFinal_days*24*3600/tStar];

opts = odeset('RelTol',1e-12, 'AbsTol',1e-12, ...
              'Events', @(t,X) earth_keepout_event_cart(t, X, rSafe));

%% Integrate
[t, X, te, Xe, ie] = ode45(@(t,X) cart_lyap_dyn_augpen(t, X, K, P, mu_nd, uMax, ...
                                    rMin, epsPenalty, wPenalty, kPenalty), ...
                           tSpan, X0, opts);

%% Recover histories
rC_nd = X(:,1:3);
vC_nd = X(:,4:6);

rT_nd_hist = X(:,7:9);
vT_nd_hist = X(:,10:12);

rC_km = rC_nd*lStar;
vC_km_s = vC_nd*vStar;

rT_km = rT_nd_hist*lStar;
vT_km_s = vT_nd_hist*vStar;

dr_km = (rC_nd - rT_nd_hist)*lStar;
dv_km_s = (vC_nd - vT_nd_hist)*vStar;

t_days = t*tStar/(24*3600);
t_sec  = t*tStar;

%% Post-processing
num_steps = length(t);

u_hist = zeros(num_steps,3);
u_norm = zeros(num_steps,1);

rmag_hist = zeros(num_steps,1);
alt_hist  = zeros(num_steps,1);

Vtrack_hist = zeros(num_steps,1);
Ppen_hist   = zeros(num_steps,1);
Vaug_hist   = zeros(num_steps,1);
g_hist      = zeros(num_steps,1);

for i = 1:num_steps

    r_i = X(i,1:3)';
    rmag_hist(i) = norm(r_i)*lStar;
    alt_hist(i) = rmag_hist(i) - rEarth;

    [~, u_out, aux] = cart_lyap_dyn_augpen_plot_helper(t(i), X(i,:)', K, P, ...
        mu_nd, uMax, rMin, epsPenalty, wPenalty, kPenalty);

    u_hist(i,:) = u_out'*aStar;
    u_norm(i) = norm(u_hist(i,:));

    Vtrack_hist(i) = aux.Vtrack;
    Ppen_hist(i)   = aux.Ppen;
    Vaug_hist(i)   = aux.Vaug;
    g_hist(i)      = aux.g;
end

%% Delta-V
DeltaV_km_s = trapz(t_sec, u_norm);
DeltaV_m_s = 1000*DeltaV_km_s;

fprintf('\nTotal Delta-V = %.6f m/s\n', DeltaV_m_s);
fprintf('Final relative position = %.6f km\n', norm(dr_km(end,:)));
fprintf('Final relative speed    = %.9f km/s\n', norm(dv_km_s(end,:)));
fprintf('Minimum altitude        = %.6f km\n', min(alt_hist));

if ~isempty(ie)
    fprintf('Terminated by event ID = %d\n', ie(end));
    if ie(end) == 1
        fprintf('Reason: Earth keep-out boundary crossed.\n');
    end
else
    fprintf('Integration ended at final time with no terminal event.\n');
end

%% Orbit plot
figure('Color','w');
hold on; grid on; axis equal; view(3);

[sx, sy, sz] = sphere(50);
topo = load('topo.mat');

surf(lStar*sx, lStar*sy, lStar*sz, ...
    'FaceColor','texturemap', ...
    'CData',topo.topo, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');

surf(rSafe*lStar*sx, rSafe*lStar*sy, rSafe*lStar*sz, ...
    'FaceAlpha',0.08, ...
    'EdgeColor','none', ...
    'FaceColor',[1 0 0], ...
    'HandleVisibility','off');

surf(rOn*lStar*sx, rOn*lStar*sy, rOn*lStar*sz, ...
    'FaceAlpha',0.04, ...
    'EdgeColor','none', ...
    'FaceColor',[0 0 0], ...
    'HandleVisibility','off');

% Initial orbit
ta_range = linspace(0, 360, 300);
r_initial = zeros(length(ta_range),3);

for j = 1:length(ta_range)
    [x, y, z] = kep2cart(a0, e0, inc0, argp0, raan0, ta_range(j), mu);
    r_initial(j,:) = [x, y, z];
end

plot3(r_initial(:,1), r_initial(:,2), r_initial(:,3), 'k:', ...
      'LineWidth',1, 'DisplayName','Initial Orbit');

% Target orbit
r_final = zeros(length(ta_range),3);

for j = 1:length(ta_range)
    [x, y, z] = kep2cart(aT, eT, incT, argpT, raanT, ta_range(j), mu);
    r_final(j,:) = [x, y, z];
end

plot3(r_final(:,1), r_final(:,2), r_final(:,3), 'b--', ...
      'LineWidth',1.2, 'DisplayName','Target Orbit');

plot3(rC_km(:,1), rC_km(:,2), rC_km(:,3), 'r', ...
      'LineWidth',1.5, 'DisplayName','Controlled Trajectory');

plot3(rT_km(:,1), rT_km(:,2), rT_km(:,3), 'b', ...
      'LineWidth',1.0, 'DisplayName','Target Trajectory');

plot3(rC_km(1,1), rC_km(1,2), rC_km(1,3), 'ro', ...
      'MarkerFaceColor','r', 'DisplayName','Chaser Start');

plot3(rT_km(1,1), rT_km(1,2), rT_km(1,3), 'bo', ...
      'MarkerFaceColor','b', 'DisplayName','Target Start');

plot3(rC_km(end,1), rC_km(end,2), rC_km(end,3), 'rs', ...
      'MarkerFaceColor','r', 'DisplayName','Chaser End');

plot3(rT_km(end,1), rT_km(end,2), rT_km(end,3), 'bs', ...
      'MarkerFaceColor','b', 'DisplayName','Target End');

xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Cartesian Lyapunov Control with Exponential Augmented Penalty');
legend('show','Location','best');

%% Cartesian state deviation
figure('Color','w');
plot(t_days, dr_km, 'LineWidth',1.5); hold on;
plot(t_days, dv_km_s, '--', 'LineWidth',1.5);
grid on;
xlabel('t [days]');
ylabel('Deviation');
title('Cartesian State Deviation');
legend('\Deltax [km]','\Deltay [km]','\Deltaz [km]', ...
       '\Deltav_x [km/s]','\Deltav_y [km/s]','\Deltav_z [km/s]', ...
       'Location','best');

%% Relative state norms
figure('Color','w');
plot(t_days, vecnorm(dr_km,2,2), 'LineWidth',1.5); hold on;
plot(t_days, vecnorm(dv_km_s,2,2), 'LineWidth',1.5);
grid on;
xlabel('t [days]');
ylabel('Norm');
title('Relative State Norms');
legend('||\Deltar|| [km]', '||\Deltav|| [km/s]', 'Location','best');

%% Control history
figure('Color','w');
plot(t_days, 1000*u_hist, 'LineWidth',1.2); hold on;
plot(t_days, 1000*u_norm, 'k', 'LineWidth',1.5);
grid on;
xlabel('t [days]');
ylabel('u [m/s^2]');
title(sprintf('Control History — \\DeltaV = %.3f m/s', DeltaV_m_s));
legend('u_x','u_y','u_z','||u||_2','Location','best');

%% Altitude history
figure('Color','w');
plot(t_days, alt_hist, 'LineWidth',1.5); hold on;
yline(hSafe_km, 'r--', 'LineWidth',1.5, 'DisplayName','Safe Altitude');
yline(hOn_km, 'k--', 'LineWidth',1.2, 'DisplayName','Penalty On Altitude');
grid on;
xlabel('t [days]');
ylabel('Altitude [km]');
title('Altitude History');
legend('Trajectory Altitude','Safe Altitude','Penalty On','Location','best');

%% Penalty history
figure('Color','w');
plot(t_days, g_hist, 'LineWidth',1.5); hold on;
yline(0, 'r--', 'LineWidth',1.2);
grid on;
xlabel('t [days]');
ylabel('g = r_{min} - ||r||');
title('Constraint Function History');
legend('g','constraint boundary','Location','best');

figure('Color','w');
semilogy(t_days, Ppen_hist, 'LineWidth',1.5);
grid on;
xlabel('t [days]');
ylabel('P(g)');
title('Exponential Penalty History');

figure('Color','w');
plot(t_days, Vtrack_hist, 'LineWidth',1.5); hold on;
plot(t_days, Vaug_hist, 'LineWidth',1.5);
grid on;
xlabel('t [days]');
ylabel('Lyapunov value');
title('Tracking Lyapunov vs Augmented Lyapunov');
legend('V_{track}','V_{aug}','Location','best');

%% ========================= FUNCTIONS =========================

function Xdot = cart_lyap_dyn_augpen(~, X, K, P, mu, uMax, ...
                                     rMin, epsPenalty, wPenalty, kPenalty)

    rC = X(1:3);
    vC = X(4:6);
    rT = X(7:9);
    vT = X(10:12);

    dr = rC - rT;
    dv = vC - vT;

    gC = -mu*rC/norm(rC)^3;
    gT = -mu*rT/norm(rT)^3;
    delta_g = gC - gT;

    % Base tracking Lyapunov function:
    % V = 1/2 dr'Kdr + 1/2 dv'dv
    Vtrack = 0.5*dr'*K*dr + 0.5*dv'*dv;

    % Minimum-radius constraint:
    % g <= 0 is safe
    rmag = norm(rC);
    g = rMin - rmag;

    % Exponential penalty:
    % P(g) = exp(k*(g + eps))
    expArg = kPenalty*(g + epsPenalty);

    % Clamp exponential argument to prevent overflow
    expArg = min(max(expArg, -50), 50);

    Ppen = exp(expArg);
    dPdg = kPenalty*Ppen;

    % dg/dr = -rhat
    dgdr = -rC/rmag;

    % Augmented Lyapunov:
    % Vaug = Vtrack + w*Vtrack*Ppen
    %
    % dVaug/dr =
    % (1 + w*Ppen)*K*dr + w*Vtrack*dPdg*dgdr
    gradV_r = (1 + wPenalty*Ppen)*(K*dr) ...
              + wPenalty*Vtrack*dPdg*dgdr;

    % Control derived from augmented Lyapunov gradient
    u = -gradV_r - P*dv - delta_g;

    % Saturate commanded control
    if norm(u) > uMax
        u = uMax*u/norm(u);
    end

    rCdot = vC;
    vCdot = gC + u;

    rTdot = vT;
    vTdot = gT;

    Xdot = [rCdot; vCdot; rTdot; vTdot];
end

function [Xdot, u_dimless, aux] = cart_lyap_dyn_augpen_plot_helper(~, X, K, P, mu, ...
                                                              uMax, rMin, ...
                                                              epsPenalty, ...
                                                              wPenalty, kPenalty)

    rC = X(1:3);
    vC = X(4:6);
    rT = X(7:9);
    vT = X(10:12);

    dr = rC - rT;
    dv = vC - vT;

    gC = -mu*rC/norm(rC)^3;
    gT = -mu*rT/norm(rT)^3;
    delta_g = gC - gT;

    Vtrack = 0.5*dr'*K*dr + 0.5*dv'*dv;

    rmag = norm(rC);
    g = rMin - rmag;

    expArg = kPenalty*(g + epsPenalty);
    expArg = min(max(expArg, -50), 50);

    Ppen = exp(expArg);
    dPdg = kPenalty*Ppen;

    dgdr = -rC/rmag;

    gradV_r = (1 + wPenalty*Ppen)*(K*dr) ...
              + wPenalty*Vtrack*dPdg*dgdr;

    u = -gradV_r - P*dv - delta_g;

    if norm(u) > uMax
        u = uMax*u/norm(u);
    end

    u_dimless = u;

    rCdot = vC;
    vCdot = gC + u;

    rTdot = vT;
    vTdot = gT;

    Xdot = [rCdot; vCdot; rTdot; vTdot];

    aux.Vtrack = Vtrack;
    aux.Ppen = Ppen;
    aux.Vaug = Vtrack + wPenalty*Vtrack*Ppen;
    aux.g = g;
end

function [value, isterminal, direction] = earth_keepout_event_cart(~, X, rSafe)

    rC = X(1:3);
    rmag = norm(rC);

    value = rmag - rSafe;
    isterminal = 1;
    direction = -1;
end