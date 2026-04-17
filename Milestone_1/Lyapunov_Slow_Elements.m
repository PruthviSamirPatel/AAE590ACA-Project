clc; clear; close all

%% Load environment
load("Earth_params.mat")
rEarth = Earth.radius;          % km
mu     = Earth.mu;              % km^3/s^2

%% Load target and parking orbit
load("Orbital_Slot.mat")
load("Parking_Orbit.mat")

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
uMax_km_s2 = 1e-5; % 0.1 m/s^2 in km/s^2
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
inc0  = deg2rad(Parking.inc);         
raan0 = deg2rad(Parking.raan);          
argp0 = deg2rad(Parking.argp);          
M0    = deg2rad(Parking.M);          

%% Target orbit
aT    = Target.a;               
eT    = Target.ecc;
incT  = deg2rad(Target.inc);        
raanT = deg2rad(Target.raan);        
argpT = deg2rad(Target.argp);        

%% Initial and target slow states (nondimensional)
xslow0_nd = [a0/lStar; e0; inc0; raan0; argp0];
xslowT_nd = [aT/lStar; eT; incT; raanT; argpT];

% Full state: [a_nd; e; i; raan; argp; M]
X0 = [xslow0_nd; M0];

%% Gains
K = eye(5);
P = eye(5);

%% Time span: at least 10 days
tFinal_days = 10;
tSpan = [0, tFinal_days*24*3600/tStar];

opts = odeset('RelTol',1e-12,'AbsTol',1e-12, ...
              'Events', @(t,X) earth_keepout_event(t, X, rSafe));

%% Integrate
[t, X] = ode45(@(t,X) kep_lyap_dyn(t, X, xslowT_nd, K, P, mu_nd, uMax, ...
                                   rSafe, rOn, kBarrier, epsBar), ...
               tSpan, X0, opts);

%% Recover histories
a_nd  = X(:,1);
ecc   = X(:,2);
inc   = X(:,3);
raan  = X(:,4);
argp  = X(:,5);
M     = X(:,6);

a_km = a_nd * lStar;

%% Post-Processing for Figures
num_steps = length(t);
r_eci = zeros(num_steps, 3);
u_hist = zeros(num_steps, 3);
u_norm = zeros(num_steps, 1);
rmag_hist = zeros(num_steps,1);
alt_hist = zeros(num_steps,1);

for i = 1:num_steps
    % current elements
    a_i = X(i,1) * lStar;
    e_i = X(i,2);
    i_i = rad2deg(X(i,3));
    O_i = rad2deg(X(i,4));
    w_i = rad2deg(X(i,5));
    M_i = rad2deg(X(i,6));
    
    % Get True Anomaly for conversion
    E_i = kepler(M_i, e_i);
    ta_i = 2*atan2(sqrt(1+e_i)*sind(E_i/2), sqrt(1-e_i)*cosd(E_i/2));
    
    % Radius magnitude
    p_i = a_i*(1-e_i^2);
    rmag_hist(i) = p_i/(1 + e_i*cos(ta_i));
    alt_hist(i) = rmag_hist(i) - rEarth;
    
    % Convert to Cartesian for Figure 22
    [x, y, z, ~, ~, ~] = kep2cart(a_i, e_i, i_i, w_i, O_i, rad2deg(ta_i), mu);
    r_eci(i,:) = [x, y, z];
    
    % Recalculate Control
    [~, u_out] = kep_lyap_dyn_plot_helper(t(i), X(i,:)', xslowT_nd, K, P, mu_nd, ...
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
ta_range = linspace(0, 2*pi, 200);
r_initial = zeros(length(ta_range), 3);
for j = 1:length(ta_range)
    [x, y, z] = kep2cart(a0, e0, rad2deg(inc0), rad2deg(argp0), ...
                         rad2deg(raan0), rad2deg(ta_range(j)), mu);
    r_initial(j,:) = [x, y, z];
end
plot3(r_initial(:,1), r_initial(:,2), r_initial(:,3), 'k:', ...
      'LineWidth', 1, 'DisplayName', 'Initial Orbit');

% Target orbit
r_final = zeros(length(ta_range), 3);
for j = 1:length(ta_range)
    [x, y, z] = kep2cart(aT, eT, rad2deg(incT), rad2deg(argpT), ...
                         rad2deg(raanT), rad2deg(ta_range(j)), mu);
    r_final(j,:) = [x, y, z];
end
plot3(r_final(:,1), r_final(:,2), r_final(:,3), 'b--', ...
      'LineWidth', 1.2, 'DisplayName', 'Target Orbit');

% Transfer trajectory
plot3(r_eci(:,1), r_eci(:,2), r_eci(:,3), 'r', ...
      'LineWidth', 1.5, 'DisplayName', 'Transfer Trajectory');

xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbit Transfer with Lyapunov Control + Earth Avoidance Barrier');
legend('show', 'Location', 'northeast');

%% Slow variable deviation
figure('Color','w');

da = X(:,1) - xslowT_nd(1);
de = X(:,2) - xslowT_nd(2);
di = zeros(length(X),1); 
dO = zeros(length(X),1); 
dw = zeros(length(X),1);

for i = 1:length(X)
    di(i) = angdiff(xslowT_nd(3), X(i,3));
    dO(i) = angdiff(xslowT_nd(4), X(i,4));
    dw(i) = angdiff(xslowT_nd(5), X(i,5));
end

plot(t_days, [da, de, di, dO, dw], 'LineWidth', 1.5);
grid on; xlabel('t [days]'); ylabel('Deviation [nonDim]');
title('Orbital Elements Deviation Over Time');
legend('sma','ecc','i','raan','argp','Location','best');

%% Control history
figure('Color','w');
plot(t_days, 1000*u_hist, 'LineWidth', 1.2); hold on;
plot(t_days, 1000*u_norm, 'k', 'LineWidth', 1.5);
grid on; xlabel('t [days]'); ylabel('u [m/s^2]');
title('Lyapunov Control History');
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

%% ========================= FUNCTIONS =========================

function [Xdot] = kep_lyap_dyn(~, X, xslowT_nd, K, P, mu, uMax, ...
                               rSafe, rOn, kBarrier, epsBar)

    % Unpack state
    a    = X(1); % nondim
    ecc  = X(2);
    inc  = X(3);
    raan = X(4);
    argp = X(5);
    M    = X(6); % rad

    if (ecc < 0 || ecc >= 1)
        error("Eccentricity out of valid elliptic range.")
    end

    % Prevent singularities from element set
    ecc_eff = max(ecc, 1e-6);
    inc_eff = inc;
    if abs(sin(inc_eff)) < 1e-6
        inc_eff = sign(inc_eff + 1e-12)*1e-6;
    end
    if abs(tan(inc_eff)) < 1e-6
        inc_eff = sign(inc_eff + 1e-12)*1e-6;
    end

    % True anomaly
    E  = kepler(rad2deg(M), ecc); % deg
    ta = 2*atan2(sqrt(1+ecc)*sind(E/2), sqrt(1-ecc)*cosd(E/2)); % rad

    % Orbital quantities
    n = sqrt(mu/a^3);
    b = a*sqrt(1-ecc^2);
    p = a*(1-ecc^2);
    h = sqrt(mu*p);
    r = p/(1 + ecc*cos(ta));

    % Unperturbed motion
    f0 = [0;0;0;0;0;n];
    
    % Gauss variational matrix
    B = 1/h * [ ...
        2*a^2*ecc*sin(ta), 2*a^2*p/r, 0;                                      % a
        p*sin(ta), (p+r)*cos(ta) + r*ecc, 0;                                  % e
        0, 0, r*cos(ta+argp);                                                 % i
        0, 0, r*sin(ta+argp)/sin(inc_eff);                                    % RAAN
        -p*cos(ta)/ecc_eff, (p+r)*sin(ta)/ecc_eff, -r*sin(ta+argp)/tan(inc_eff); % argp
        (b*p*cos(ta))/(a*ecc_eff) - (2*b*r)/a, -b*(p+r)*sin(ta)/(a*ecc_eff), 0   % M
        ];

    B_slow = B(1:5, :);

    % Slow-state error
    xslow = [a; ecc; inc; raan; argp];
    dx_slow = zeros(5,1);
    dx_slow(1) = xslow(1) - xslowT_nd(1);
    dx_slow(2) = xslow(2) - xslowT_nd(2);
    dx_slow(3) = angdiff(xslowT_nd(3),  xslow(3));
    dx_slow(4) = angdiff(xslowT_nd(4),  xslow(4));
    dx_slow(5) = angdiff(xslowT_nd(5),  xslow(5));

    % Nominal Lyapunov controller
    H = B_slow' * K' * K * B_slow;
    g = B_slow' * K' * P * dx_slow;
    uLyap = -(H \ g);

    % Barrier / penalty control
    d   = r - rSafe;   % clearance above safety radius
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        % Emergency outward radial push
        uBarrier(1) = kBarrier / epsBar;
    elseif d < dOn
        % Smooth repulsive barrier that is zero at r = rOn
        uBarrier(1) = kBarrier * (1/(d + epsBar) - 1/dOn);
    end

    % Total control
    u = uLyap + uBarrier;

    % Saturate
    if norm(u) > uMax
        u = uMax * u / norm(u);
    end

    Xdot = f0 + B * u;
end

%% Plot helper
function [Xdot, u_dimless] = kep_lyap_dyn_plot_helper(~, X, xslowT_nd, K, P, ...
                                                      mu, uMax, rSafe, rOn, ...
                                                      kBarrier, epsBar)
    % Unpack state
    a    = X(1);
    ecc  = X(2);
    inc  = X(3);
    raan = X(4);
    argp = X(5);
    M    = X(6);

    if (ecc < 0 || ecc >= 1)
        error("Eccentricity out of valid elliptic range.")
    end

    % Prevent singularities
    ecc_eff = max(ecc, 1e-6);
    inc_eff = inc;
    if abs(sin(inc_eff)) < 1e-6
        inc_eff = sign(inc_eff + 1e-12)*1e-6;
    end
    if abs(tan(inc_eff)) < 1e-6
        inc_eff = sign(inc_eff + 1e-12)*1e-6;
    end

    % True anomaly
    E  = kepler(rad2deg(M), ecc);
    ta = 2*atan2(sqrt(1+ecc)*sind(E/2), sqrt(1-ecc)*cosd(E/2));

    % Orbital quantities
    n = sqrt(mu/a^3);
    b = a*sqrt(1-ecc^2);
    p = a*(1-ecc^2);
    h = sqrt(mu*p);
    r = p/(1 + ecc*cos(ta));

    % Unperturbed motion
    f0 = [0;0;0;0;0;n];

    % Gauss matrix
    B = 1/h * [ ...
        2*a^2*ecc*sin(ta), 2*a^2*p/r, 0;
        p*sin(ta), (p+r)*cos(ta) + r*ecc, 0;
        0, 0, r*cos(ta+argp);
        0, 0, r*sin(ta+argp)/sin(inc_eff);
        -p*cos(ta)/ecc_eff, (p+r)*sin(ta)/ecc_eff, -r*sin(ta+argp)/tan(inc_eff);
        (b*p*cos(ta))/(a*ecc_eff) - (2*b*r)/a, -b*(p+r)*sin(ta)/(a*ecc_eff), 0
        ];

    B_slow = B(1:5, :);

    % Slow-state error
    xslow = [a; ecc; inc; raan; argp];
    dx_slow = zeros(5,1);
    dx_slow(1) = xslow(1) - xslowT_nd(1);
    dx_slow(2) = xslow(2) - xslowT_nd(2);
    dx_slow(3) = angdiff(xslowT_nd(3),  xslow(3));
    dx_slow(4) = angdiff(xslowT_nd(4),  xslow(4));
    dx_slow(5) = angdiff(xslowT_nd(5),  xslow(5));

    % Nominal Lyapunov control
    H = B_slow' * K' * K * B_slow;
    g = B_slow' * K' * P * dx_slow;
    uLyap = -(H \ g);

    % Barrier control
    d   = r - rSafe;
    dOn = rOn - rSafe;

    uBarrier = [0;0;0];

    if d <= 0
        uBarrier(1) = kBarrier / epsBar;
    elseif d < dOn
        uBarrier(1) = kBarrier * (1/(d + epsBar) - 1/dOn);
    end

    % Total control
    u = uLyap + uBarrier;

    % Saturate
    if norm(u) > uMax
        u = uMax * u / norm(u);
    end

    u_dimless = u;
    Xdot = f0 + B*u;
end

%% Event function: terminate if the safe radius is crossed
function [value, isterminal, direction] = earth_keepout_event(~, X, rSafe)
    a   = X(1);
    ecc = X(2);
    M   = X(6);

    % Handle bad values robustly
    if ecc < 0 || ecc >= 1
        value = -1;
        isterminal = 1;
        direction = -1;
        return
    end

    E  = kepler(rad2deg(M), ecc);
    ta = 2*atan2(sqrt(1+ecc)*sind(E/2), sqrt(1-ecc)*cosd(E/2));

    p = a*(1-ecc^2);
    r = p/(1 + ecc*cos(ta));

    value = r - rSafe;
    isterminal = 1;
    direction = -1;
end