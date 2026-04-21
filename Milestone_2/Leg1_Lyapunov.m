clc; clear; close all

%% Constants
days_coast = 10;
t_coast = days_coast*24*60*60; % s

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
t_coast_nd = t_coast/tStar;

%% Max control
uMax_km_s2 = 1e-4; % 10 cm/s = 1e-5km/s
uMax = uMax_km_s2 / (lStar / tStar^2); % nondimensionalized

%% Earth avoidance / safety barrier settings
hSafe_km = 150;                 % minimum allowed altitude above Earth surface
hOn_km   = 300;                 % altitude where barrier starts turning on
rSafe    = (rEarth + hSafe_km) / lStar;   % nondim safe radius
rOn      = (rEarth + hOn_km)   / lStar;   % nondim activation radius
kBarrier = 5e-5;               % barrier gain
epsBar   = 1e-6;               % avoids singularity in barrier

%% Initial orbit (classical)
a0    = Parking.a;
e0    = Parking.ecc;
inc0  = deg2rad(Parking.inc);
raan0 = deg2rad(Parking.raan);
argp0 = deg2rad(Parking.argp);
M0    = deg2rad(Parking.M);

%% Intermediary orbit (classical)
Inter = get_intermediary_orbit(Earth, Parking, Target, t_coast);
aI    = Inter.a;
eI    = Inter.ecc;
incI  = deg2rad(Inter.inc);
raanI = deg2rad(Inter.raan);
argpI = deg2rad(Inter.argp);

%% Initial true anomaly from mean anomaly
E0 = kepler(rad2deg(M0), e0); % deg
ta0 = 2*atan2(sqrt(1+e0)*sind(E0/2), sqrt(1-e0)*cosd(E0/2)); % rad

%% Target orbit: choose TA = 0 as reference shape/orientation target
% For slow elements, TA is not part of the target Lyapunov state anyway.
taI = 0;

%% Convert initial and target to MEE (dimensional first)
mee0_dim = kep2mee(a0, e0, inc0, raan0, argp0, ta0);
meeT_dim = kep2mee(aI, eI, incI, raanI, argpI, taI);

%% Nondimensionalize propagated states
X0 = [mee0_dim(1)/lStar; mee0_dim(2:6)];
xslowT_nd = [meeT_dim(1)/lStar; meeT_dim(2:5)];

%% Gains
K = eye(5);
P = eye(5);

%% Time span
tFinal_days = 0.25;
tSpan = [0, tFinal_days*24*3600/tStar];

opts = odeset('RelTol',1e-12,'AbsTol',1e-12, ...
              'Events', @(t,X) earth_keepout_event_mee(t, X, rSafe));

%% Integrate
Re_nd = Earth.radius / lStar;
J2 = Earth.J2;

[t, X] = ode45(@(t,X) mee_lyap_dyn(t, X, xslowT_nd, K, P, mu_nd, uMax, ...
                                   rSafe, rOn, kBarrier, epsBar, J2, Re_nd), ...
               tSpan, X0, opts);
%% Recover histories
num_steps = length(t);

p_nd = X(:,1);
f_hist = X(:,2);
g_hist = X(:,3);
h_hist = X(:,4);
k_hist = X(:,5);
L_hist = X(:,6);

a_km      = zeros(num_steps,1);
ecc_hist  = zeros(num_steps,1);
inc_hist  = zeros(num_steps,1);
raan_hist = zeros(num_steps,1);
argp_hist = zeros(num_steps,1);
ta_hist   = zeros(num_steps,1);

r_eci = zeros(num_steps, 3);
u_hist = zeros(num_steps, 3);
u_norm = zeros(num_steps, 1);
rmag_hist = zeros(num_steps,1);
alt_hist = zeros(num_steps,1);

for i = 1:num_steps
    mee_dim_i = [X(i,1)*lStar; X(i,2:6)'];
    [a_i, e_i, inc_i, raan_i, argp_i, ta_i] = mee2kep(mee_dim_i);

    a_km(i)      = a_i;
    ecc_hist(i)  = e_i;
    inc_hist(i)  = inc_i;
    raan_hist(i) = raan_i;
    argp_hist(i) = argp_i;
    ta_hist(i)   = ta_i;

    % Radius from MEE directly
    p_i_nd = X(i,1);
    f_i = X(i,2);
    g_i = X(i,3);
    L_i = X(i,6);
    w_i = 1 + f_i*cos(L_i) + g_i*sin(L_i);
    rmag_hist(i) = (p_i_nd / w_i) * lStar;
    alt_hist(i) = rmag_hist(i) - rEarth;

    % Convert to Cartesian
    [x, y, z, ~, ~, ~] = kep2cart(a_i, e_i, rad2deg(inc_i), ...
                                  rad2deg(argp_i), rad2deg(raan_i), ...
                                  rad2deg(ta_i), mu);
    r_eci(i,:) = [x, y, z];

    % Recompute control for plotting
    [~, u_out] = mee_lyap_dyn_plot_helper(t(i), X(i,:)', xslowT_nd, K, P, mu_nd, ...
                                          uMax, rSafe, rOn, kBarrier, epsBar, J2, Re_nd);
    u_hist(i,:) = u_out' * (lStar / tStar^2);
    u_norm(i) = norm(u_hist(i,:));
end

t_days = (t * tStar) / (24*3600);

%% Orbit plot
figure('Color','w');
hold on; grid on; axis equal; view(3);

[sx, sy, sz] = sphere(50);
topo = load('topo.mat');
surf(lStar*sx, lStar*sy, lStar*sz, 'FaceColor', 'texturemap', ...
     'CData', topo.topo, 'EdgeColor', 'none', 'HandleVisibility','off');

surf(rSafe*lStar*sx, rSafe*lStar*sy, rSafe*lStar*sz, ...
    'FaceAlpha', 0.08, 'EdgeColor', 'none', 'FaceColor', [1 0 0], ...
    'HandleVisibility','off');

% Initial orbit
ta_range = linspace(0, 2*pi, 300);
r_initial = zeros(length(ta_range), 3);
for j = 1:length(ta_range)
    [x, y, z] = kep2cart(a0, e0, rad2deg(inc0), rad2deg(argp0), ...
                         rad2deg(raan0), rad2deg(ta_range(j)), mu);
    r_initial(j,:) = [x, y, z];
end
plot3(r_initial(:,1), r_initial(:,2), r_initial(:,3), 'k:', ...
      'LineWidth', 1, 'DisplayName', 'Initial Orbit');

% Inter orbit
r_final = zeros(length(ta_range), 3);
for j = 1:length(ta_range)
    [x, y, z] = kep2cart(aI, eI, rad2deg(incI), rad2deg(argpI), ...
                         rad2deg(raanI), rad2deg(ta_range(j)), mu);
    r_final(j,:) = [x, y, z];
end
plot3(r_final(:,1), r_final(:,2), r_final(:,3), 'b--', ...
      'LineWidth', 1.2, 'DisplayName', 'Inter Orbit');

% Transfer trajectory
plot3(r_eci(:,1), r_eci(:,2), r_eci(:,3), 'r', ...
      'LineWidth', 1.5, 'DisplayName', 'Transfer Trajectory');

xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbit Transfer with Lyapunov Control in MEE + Earth Avoidance Barrier');
legend('show', 'Location', 'northeast');

%% Slow variable deviation
figure('Color','w');
dp = X(:,1) - xslowT_nd(1);
df = X(:,2) - xslowT_nd(2);
dg = X(:,3) - xslowT_nd(3);
dh = X(:,4) - xslowT_nd(4);
dk = X(:,5) - xslowT_nd(5);

plot(t_days, [dp, df, dg, dh, dk], 'LineWidth', 1.5);
grid on; xlabel('t [days]'); ylabel('Deviation [nonDim]');
title('MEE Slow-State Deviation Over Time');
legend('p','f','g','h','k','Location','best');

%% Control history
t_sec = t * tStar;          % convert time to seconds
DeltaV_km_s = trapz(t_sec, u_norm);

figure('Color','w');
plot(t_days, 1000*u_hist, 'LineWidth', 1.2); hold on;
plot(t_days, 1000*u_norm, 'k', 'LineWidth', 1.5);
grid on; xlabel('t [days]'); ylabel('u [m/s^2]');
title(sprintf('Lyapunov Control History — ∆v = %.3f km/s', DeltaV_km_s));
legend('u_r','u_t','u_n','||u||_2','Location','best');

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

%% Classical element history
figure('Color','w');
subplot(3,2,1)
plot(t_days, a_km, 'LineWidth', 1.5); grid on
xlabel('t [days]'); ylabel('a [km]')

subplot(3,2,2)
plot(t_days, ecc_hist, 'LineWidth', 1.5); grid on
xlabel('t [days]'); ylabel('e')

subplot(3,2,3)
plot(t_days, rad2deg(inc_hist), 'LineWidth', 1.5); grid on
xlabel('t [days]'); ylabel('i [deg]')

subplot(3,2,4)
plot(t_days, rad2deg(raan_hist), 'LineWidth', 1.5); grid on
xlabel('t [days]'); ylabel('\Omega [deg]')

subplot(3,2,5)
plot(t_days, rad2deg(argp_hist), 'LineWidth', 1.5); grid on
xlabel('t [days]'); ylabel('\omega [deg]')

subplot(3,2,6)
plot(t_days, rad2deg(ta_hist), 'LineWidth', 1.5); grid on
xlabel('t [days]'); ylabel('\nu [deg]')

sgtitle('Recovered Classical Orbital Elements')

%% Save final orbital elements
InterAchieved.a = a_km(end);
InterAchieved.ecc = ecc_hist(end);
InterAchieved.inc = rad2deg(inc_hist(end));
InterAchieved.raan = rad2deg(raan_hist(end));
InterAchieved.argp = rad2deg(argp_hist(end));

% Final mean anomaly
ta_end = ta_hist(end);              % rad
e_end  = ecc_hist(end);

E_end = 2*atan2( sqrt(1-e_end)*sin(ta_end/2), ...
                 sqrt(1+e_end)*cos(ta_end/2) );   % rad
E_end = mod(E_end, 2*pi);

M_end = E_end - e_end*sin(E_end);   % rad
M_end = mod(M_end, 2*pi);

InterAchieved.M = rad2deg(M_end);

% Mean motion
InterAchieved.n = sqrt(Earth.mu / InterAchieved.a^3);   % rad/s
save('Intermediate_Orbit.mat', 'InterAchieved')

%% Plot helper
function [Xdot, u_dimless] = mee_lyap_dyn_plot_helper(~, X, xslowT_nd, K, P, ...
                                                      mu, uMax, rSafe, rOn, ...
                                                      kBarrier, epsBar, J2, Re_nd)

    p = X(1);
    f = X(2);
    g = X(3);
    h = X(4);
    k = X(5);
    L = X(6);

    if p <= 0 || ~isfinite(p)
        Xdot = zeros(6,1);
        u_dimless = zeros(3,1);
        return
    end

    w  = 1 + f*cos(L) + g*sin(L);
    s2 = 1 + h^2 + k^2;
    r  = p / w;

    if w <= 1e-10 || r <= 0 || ~isfinite(r)
        Xdot = zeros(6,1);
        u_dimless = zeros(3,1);
        return
    end

    sq = sqrt(p/mu);

    f0 = [0; 0; 0; 0; 0; sqrt(mu/p^3)*w^2];

    B = zeros(6,3);

    B(1,:) = [0,            2*p/w*sq,                          0];

    B(2,:) = [sq*sin(L),    sq*((w+1)*cos(L) + f)/w, ...
                         -sq*(h*sin(L) - k*cos(L))/w];

    B(3,:) = [-sq*cos(L),   sq*((w+1)*sin(L) + g)/w, ...
                          sq*(h*cos(L) + k*sin(L))/w];

    B(4,:) = [0,            0,  sq*(s2/(2*w))*cos(L)];

    B(5,:) = [0,            0,  sq*(s2/(2*w))*sin(L)];

    B(6,:) = [0,            0,  sq*(h*sin(L) - k*cos(L))/w];

    Bslow = B(1:5,:);

    xslow = X(1:5);
    dx_slow = xslow - xslowT_nd;

    H = Bslow' * K' * K * Bslow + 1e-12*eye(3);
    gvec = Bslow' * K' * P * dx_slow;
    uLyap = -(H \ gvec);

    d   = r - rSafe;
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        uBarrier(1) = kBarrier / epsBar;
    elseif d < dOn
        uBarrier(1) = kBarrier * (1/(d + epsBar) - 1/dOn);
    end

    u = uLyap + uBarrier;

    unorm = norm(u);
    if unorm > uMax
        u = uMax * u / unorm;
    end

    u_dimless = u;
    [r_eci, v_eci] = mee2rv(X, mu);
    aJ2_eci = accelJ2(r_eci, mu, J2, Re_nd);
    aJ2_rtn = eci2rtn_accel(r_eci, v_eci, aJ2_eci);
    
    Xdot = f0 + B*(u + aJ2_rtn);
end

%% Event function
function [value, isterminal, direction] = earth_keepout_event_mee(~, X, rSafe)
    p = X(1);
    f = X(2);
    g = X(3);
    L = X(6);

    w = 1 + f*cos(L) + g*sin(L);

    if p <= 0 || w <= 0
        value = -1;
        isterminal = 1;
        direction = -1;
        return
    end

    r = p / w;

    value = r - rSafe;
    isterminal = 1;
    direction = -1;
end