clc; clear; close all;
% this script solves the minimum time optimal control transfer problem from
% parking orbit to target slot. It assumes no perturbations. Particularly,
% assumes that target orbit does not get purturbed.

%% Import Earth environment
load("Earth_params.mat")
rEarth = Earth.radius; 
mu = Earth.mu;

%% Load Target Orbit
load('Orbital_Slot.mat')

%% Get nondim quantities
lStar = rEarth;
tStar = sqrt(rEarth^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;
mu_nd = 1;

%% Constants
umax_dim = 1e-3; % 1 m/s2
umax = umax_dim / aStar;
B = [zeros(3,3); eye(3)];
t0 = 0;
opts = odeset('RelTol',1e-12, 'AbsTol',1e-12);
numDays = .2;
tf = numDays*24*60*60;
tf = tf/tStar;

%% Load Parking Orbit
load('Parking_Orbit.mat')

% get true anomaly:
E = kepler(Parking.M,Parking.ecc); 
ta = 2 * atan2(sqrt(1 + Parking.ecc) * sind(E/2), ...
               sqrt(1 - Parking.ecc) * cosd(E/2));
ta = rad2deg(ta);
% get intial cartesian state:
[x0, y0, z0, vx0, vy0, vz0] = kep2cart(Parking.a, Parking.ecc, Parking.inc, ...
    Parking.argp, Parking.raan, ta, mu);

% nondimensionalize:
r0 = [x0; y0; z0]/lStar;
v0 = [vx0; vy0; vz0]/vStar;
X0 = [r0;v0];

%% Get orbit slot state (assuming no perturbations) at tf
% mean anomaly at tf:
M_f = Target.M + rad2deg(Target.n_nd * (tf - t0));
M_f = mod(M_f, 360);

% convert to true anomaly:
E = kepler(M_f,Target.ecc); 
ta = 2 * atan2(sqrt(1 + Target.ecc) * sind(E/2), ...
               sqrt(1 - Target.ecc) * cosd(E/2));
ta = rad2deg(ta);
% get target cartesian state:
[xTf, yTf, zTf, vxTf, vyTf, vzTf] = kep2cart(Target.a_nd, Target.ecc, Target.inc, ...
    Target.argp, Target.raan, ta, mu_nd);
xTargetf = [xTf; yTf; zTf; vxTf; vyTf; vzTf];

%% get initial lambda with fsolve
fsolve_opts = optimoptions('fsolve', 'OptimalityTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, 'Display', 'iter');

rhos = [1, .5, .2, .1, 0.05, 1e-2];
lambda0s = zeros(6, length(rhos));
% lambda0_guess = [0.8; 0.1; 0.2; 1.1];
lambda0_guess = [1e-3; -1e-3; 1e-3; -1e-3; 1e-3; -1e-3];
for i=1:length(rhos)
    rho = rhos(i);
    fun = @(lambda0) shooting_minFuel(lambda0, X0, t0, tf, umax, mu_nd, B, xTargetf, rho, opts); 
    [lambda0, fval, exitflag] = fsolve(fun, lambda0_guess, fsolve_opts);
    lambda0_guess = lambda0;
    lambda0s(:,i) = lambda0;
end
fprintf('\nConverged Initial Costate: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', lambda0);

%% Generate Optimal Trajectory for Plotting
N = 1000;
t = linspace(t0, tf, N);

X_init = [X0; lambda0];
sol = ode45(@(t,X) EoM_minFuel(t, X, B, mu_nd, umax, rho), [t0 tf], X_init, opts);
X_eval = deval(sol, t);

pos = X_eval(1:3, :);       % nondimensional transfer position
vel = X_eval(4:6, :);
lambda = X_eval(7:12, :);

fprintf("\n\n Position Difference at tf: %e", norm(pos(:,end) - xTargetf(1:3)))
fprintf("\n\n Velocity Difference at tf: %e", norm(vel(:,end) - xTargetf(4:6)))

% dimensional transfer position for plotting
pos_dim = pos * lStar;

%% Calculate Control History
u = zeros(3, N);
for k = 1:N
    p = -B' * lambda(:,k);
    uHat = p / norm(p);
    
    % Switching Function: S > 0 when ||p|| > 1
    S = norm(p) - 1; 
    
    Gamma = umax/2 * (1 + tanh(S/rho));
    u = uHat * Gamma;

    if norm(p) < 1e-10
        u(:,k) = zeros(3,1);
    else
        u(:,k) = u;
    end
end

%% Generate Initial Parking Orbit for Plotting
Tpark_nd = 2*pi*sqrt((Parking.a/lStar)^3 / mu_nd);
tPark = linspace(0, Tpark_nd, N);

park_state0_nd = [r0; v0];
solPark = ode45(@(t,X) twoBodyEOM_nd(t, X, mu_nd), [0 Tpark_nd], park_state0_nd, opts);
Xpark = deval(solPark, tPark);

rPark_dim = Xpark(1:3,:) * lStar;

%% Generate Target State Evolution for Plotting
target_hist = zeros(6, N);

for k = 1:N
    M_k = Target.M + rad2deg(Target.n_nd * (t(k) - t0));
    M_k = mod(M_k, 360);

    E_k = kepler(M_k, Target.ecc);
    ta_k = 2 * atan2(sqrt(1 + Target.ecc) * sind(E_k/2), ...
                     sqrt(1 - Target.ecc) * cosd(E_k/2));
    ta_k = rad2deg(ta_k);

    [xT, yT, zT, vxT, vyT, vzT] = kep2cart(Target.a, Target.ecc, Target.inc, ...
        Target.argp, Target.raan, ta_k, mu);

    target_hist(:,k) = [xT; yT; zT; vxT; vyT; vzT];
end

rTarget_dim = target_hist(1:3,:);

%% Earth sphere for plotting
[xe, ye, ze] = sphere(60);

%% 3D Trajectory Plot
figure;
surf(rEarth*xe, rEarth*ye, rEarth*ze, ...
    'FaceColor', [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;
grid on;
axis equal;
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Minimum-Time Transfer: Parking Orbit to Target Slot');

plot3(rPark_dim(1,:), rPark_dim(2,:), rPark_dim(3,:), 'k--', 'LineWidth', 1.2);
plot3(rTarget_dim(1,:), rTarget_dim(2,:), rTarget_dim(3,:), 'r-', 'LineWidth', 1.5);
plot3(pos_dim(1,:), pos_dim(2,:), pos_dim(3,:), 'b-', 'LineWidth', 2);

plot3(rPark_dim(1,1), rPark_dim(2,1), rPark_dim(3,1), 'ko', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 6);
plot3(rTarget_dim(1,1), rTarget_dim(2,1), rTarget_dim(3,1), 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 6);
plot3(pos_dim(1,end), pos_dim(2,end), pos_dim(3,end), 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 6);

legend('Earth', 'Initial parking orbit', 'Target state evolution', ...
       'Optimal transfer arc', 'Initial spacecraft state', ...
       'Initial target state', 'Final transfer state', ...
       'Location', 'best');

view(3);

%% Optional: Position vs Time
figure;
plot(t*tStar/3600, vecnorm(pos_dim,2,1), 'b', 'LineWidth', 1.5); hold on;
plot(t*tStar/3600, vecnorm(rTarget_dim,2,1), 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time [hr]');
ylabel('Radius [km]');
title('Radius History');
legend('Transfer trajectory', 'Target evolution', 'Location', 'best');

%% Optional: Control magnitude vs Time
figure;
plot(t*tStar/3600, vecnorm(u,2,1), 'LineWidth', 1.5);
grid on;
xlabel('Time [hr]');
ylabel('||u|| [nondim]');
title('Control Magnitude History');

%% Fuel cost
dt = t(2) - t(1);
J = sum(vecnorm(u,2,1)) * dt

%% FUNCTION HELPER
function Xdot = twoBodyEOM_nd(~, X, mu_nd)
    r = X(1:3);
    v = X(4:6);

    r_norm = norm(r);

    a = -mu_nd * r / r_norm^3;

    Xdot = [v; a];
end




