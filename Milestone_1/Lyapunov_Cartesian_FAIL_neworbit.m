clc; clear; close all

%% Load environment
load("Earth_params.mat")
rEarth = Earth.radius;          % km
mu     = Earth.mu;              % km^3/s^2

%% Load target and parking orbit
load("Orbital_Slot_new.mat")
load("Parking_Orbit_new.mat")

%% Nondimensionalization
r0 = rEarth;
lStar = rEarth;
tStar = sqrt(rEarth^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;
mu_nd = 1;

%% Nondimensional scales
lStar = r0;
tStar = sqrt(lStar^3 / mu);
aStar = lStar / tStar^2;     % km/s^2
uMax_km_s2 = 1e-5;           % 0.1 m/s^2 in km/s^2
uMax = uMax_km_s2 / (lStar / tStar^2); % Non-dimensionalized

%% Earth avoidance / safety barrier settings
hSafe_km = 150;                 % minimum allowed altitude above Earth surface
hOn_km   = 300;                 % altitude where barrier starts turning on
rSafe    = (rEarth + hSafe_km) / lStar;   % nondim safe radius
rOn      = (rEarth + hOn_km)   / lStar;   % nondim activation radius
kBarrier = 5e-5;               % barrier gain (tune this)
epsBar   = 1e-6;               % avoids singularity in barrier

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

%% Initial and target Cartesian states (nondimensional)
[x, y, z, vx, vy, vz] = kep2cart(a0, e0, inc0, argp0, raan0, M0, mu);
rC0_nd = [x; y; z] / lStar;
vC0_nd = [vx; vy; vz] / vStar;

[x, y, z, vx, vy, vz] = kep2cart(aT, eT, incT, argpT, raanT, MT, mu);
rT_nd = [x; y; z] / lStar;
vT_nd = [vx; vy; vz] / vStar;

% Full state: [rC; vC; rT; vT]
X0 = [rC0_nd; vC0_nd; rT_nd; vT_nd];

%% Gains
K = eye(3);
P = eye(3);

%% Time span: at least 10 days
tFinal_days = 10;
tSpan = [0, tFinal_days*24*3600/tStar];

opts = odeset('RelTol',1e-12,'AbsTol',1e-12, ...
              'Events', @(t,X) earth_keepout_event_cart(t, X, rSafe));

%% Integrate
[t, X] = ode45(@(t,X) cart_lyap_dyn(t, X, K, P, mu_nd, uMax, ...
                                    rSafe, rOn, kBarrier, epsBar), ...
               tSpan, X0, opts);

%% Recover histories
rC_nd = X(:,1:3);
vC_nd = X(:,4:6);
rT_nd_hist = X(:,7:9);
vT_nd_hist = X(:,10:12);

rC_km = rC_nd * lStar;
rT_km = rT_nd_hist * lStar;

dr_km = (rC_nd - rT_nd_hist) * lStar;
dv_km_s = (vC_nd - vT_nd_hist) * vStar;

%% Post-Processing for Figures
num_steps = length(t);
u_hist = zeros(num_steps, 3);
u_norm = zeros(num_steps, 1);
rmag_hist = zeros(num_steps,1);
alt_hist = zeros(num_steps,1);

for i = 1:num_steps
    % Current controlled spacecraft state
    r_i = X(i,1:3)';
    v_i = X(i,4:6)';
    rT_i = X(i,7:9)';
    vT_i = X(i,10:12)';
    
    % Radius magnitude
    rmag_hist(i) = norm(r_i) * lStar;
    alt_hist(i) = rmag_hist(i) - rEarth;
    
    % Recalculate Control
    [~, u_out] = cart_lyap_dyn_plot_helper(t(i), X(i,:)', K, P, mu_nd, ...
                                           uMax, rSafe, rOn, kBarrier, epsBar);
    u_hist(i,:) = u_out' * (lStar / tStar^2);
    u_norm(i) = norm(u_hist(i,:));
end

t_days = (t * tStar) / (24*3600);

%% Orbit plot
figure('Color','w');
hold on; grid on; axis equal; view(3);

% Plot Earth
[sx, sy, sz] = sphere(50);
topo = load('topo.mat');
surf(lStar*sx, lStar*sy, lStar*sz, 'FaceColor', 'texturemap', ...
     'CData', topo.topo, 'EdgeColor', 'none', 'HandleVisibility','off');

% Plot safety sphere
surf(rSafe*lStar*sx, rSafe*lStar*sy, rSafe*lStar*sz, ...
    'FaceAlpha', 0.08, 'EdgeColor', 'none', 'FaceColor', [1 0 0], ...
    'HandleVisibility','off');

% Initial orbit
ta_range = linspace(0, 360, 200);
r_initial = zeros(length(ta_range), 3);
for j = 1:length(ta_range)
    [x, y, z] = kep2cart(a0, e0, inc0, argp0, raan0, ta_range(j), mu);
    r_initial(j,:) = [x, y, z];
end
plot3(r_initial(:,1), r_initial(:,2), r_initial(:,3), 'k:', ...
      'LineWidth', 1, 'DisplayName', 'Initial Orbit');

% Target orbit
r_final = zeros(length(ta_range), 3);
for j = 1:length(ta_range)
    [x, y, z] = kep2cart(aT, eT, incT, argpT, raanT, ta_range(j), mu);
    r_final(j,:) = [x, y, z];
end
plot3(r_final(:,1), r_final(:,2), r_final(:,3), 'b--', ...
      'LineWidth', 1.2, 'DisplayName', 'Target Orbit');

% Transfer trajectory
plot3(rC_km(:,1), rC_km(:,2), rC_km(:,3), 'r', ...
      'LineWidth', 1.5, 'DisplayName', 'Transfer Trajectory');

xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbit Transfer with Cartesian Lyapunov Control + Earth Avoidance Barrier');
legend('show', 'Location', 'northeast');

%% Cartesian state deviation
figure('Color','w');

plot(t_days, [dr_km, dv_km_s], 'LineWidth', 1.5);
grid on; xlabel('t [days]'); ylabel('Deviation');
title('Cartesian State Deviation Over Time');
legend('\Deltax [km]','\Deltay [km]','\Deltaz [km]', ...
       '\Delta v_x [km/s]','\Delta v_y [km/s]','\Delta v_z [km/s]', ...
       'Location','best');

%% Control history
figure('Color','w');
plot(t_days, 1000*u_hist, 'LineWidth', 1.2); hold on;
plot(t_days, 1000*u_norm, 'k', 'LineWidth', 1.5);
grid on; xlabel('t [days]'); ylabel('u [m/s^2]');
title('Lyapunov Control History');
legend('u_x','u_y','u_z','||u||_2','Location','best');

%% Altitude history
figure('Color','w');
plot(t_days, alt_hist, 'LineWidth', 1.5); hold on;
yline(hSafe_km, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Safe Altitude');
yline(hOn_km, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Barrier On Altitude');
grid on;
xlabel('t [days]');
ylabel('Altitude [km]');
title('Altitude History');
legend('Trajectory Altitude','Safe Altitude','Barrier On','Location','best');

%% ========================= FUNCTIONS =========================

function [Xdot] = cart_lyap_dyn(~, X, K, P, mu, uMax, ...
                                rSafe, rOn, kBarrier, epsBar)

    % Unpack state
    rC = X(1:3);
    vC = X(4:6);
    rT = X(7:9);
    vT = X(10:12);

    % Relative state
    dr = rC - rT;
    dv = vC - vT;

    % Two-body gravity
    gC = -mu * rC / norm(rC)^3;
    gT = -mu * rT / norm(rT)^3;
    delta_g = gC - gT;

    % Nominal Lyapunov controller
    uLyap = -K*dr - P*dv - delta_g;

    % Barrier / penalty control
    rmag = norm(rC);
    d   = rmag - rSafe;
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        % Emergency outward radial push
        uBarrier = (kBarrier / epsBar) * (rC / norm(rC));
    elseif d < dOn
        % Smooth repulsive barrier that is zero at r = rOn
        uBarrier = kBarrier * (1/(d + epsBar) - 1/dOn) * (rC / norm(rC));
    end

    % Total control
    u = uLyap + uBarrier;

    % Saturate
    if norm(u) > uMax
        u = uMax * u / norm(u);
    end

    % Equations of motion
    rCdot = vC;
    vCdot = gC + u;

    rTdot = vT;
    vTdot = gT;

    Xdot = [rCdot; vCdot; rTdot; vTdot];
end

%% Plot helper
function [Xdot, u_dimless] = cart_lyap_dyn_plot_helper(~, X, K, P, mu, ...
                                                       uMax, rSafe, rOn, ...
                                                       kBarrier, epsBar)

    % Unpack state
    rC = X(1:3);
    vC = X(4:6);
    rT = X(7:9);
    vT = X(10:12);

    % Relative state
    dr = rC - rT;
    dv = vC - vT;

    % Two-body gravity
    gC = -mu * rC / norm(rC)^3;
    gT = -mu * rT / norm(rT)^3;
    delta_g = gC - gT;

    % Nominal Lyapunov controller
    uLyap = -K*dr - P*dv - delta_g;

    % Barrier control
    rmag = norm(rC);
    d   = rmag - rSafe;
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        uBarrier = (kBarrier / epsBar) * (rC / norm(rC));
    elseif d < dOn
        uBarrier = kBarrier * (1/(d + epsBar) - 1/dOn) * (rC / norm(rC));
    end

    % Total control
    u = uLyap + uBarrier;

    % Saturate
    if norm(u) > uMax
        u = uMax * u / norm(u);
    end

    u_dimless = u;

    % Equations of motion
    rCdot = vC;
    vCdot = gC + u;

    rTdot = vT;
    vTdot = gT;

    Xdot = [rCdot; vCdot; rTdot; vTdot];
end

%% Event function: terminate if the safe radius is crossed
function [value, isterminal, direction] = earth_keepout_event_cart(~, X, rSafe)
    rC = X(1:3);
    rmag = norm(rC);

    value = rmag - rSafe;
    isterminal = 1;
    direction = -1;
end