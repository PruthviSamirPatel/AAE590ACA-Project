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

%% Solve Two-Point Boundary Value Problem (TPBVP)
% Optimization vector Z = [lambda1; lambda2; lambda3; lambda4; lambda5; lambda6; tf]
% Z_guess = [1e-3; -1e-3; 1e-3; -1e-3; 1e-3; -1e-3; 6.5]; % this works well
% to converge as an inital guess. 
Z_guess = [4.7758, 5.6636, 3.3864, -2.9888, 5.5177, 7.4891, 5.9236]';

fsolve_opts = optimoptions('fsolve', 'OptimalityTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, 'Display', 'iter');

[Z_sol, fval, exitflag] = fsolve(@(Z) shooting_minTime(Z, X0, t0, umax, mu_nd, B, Target, opts), Z_guess, fsolve_opts);

lambda0 = Z_sol(1:6);
tf = Z_sol(7);

fprintf('\nConverged Initial Costate: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', lambda0);
fprintf('Converged Final Time (tf): %.4f\n', tf);

%% Generate Optimal Trajectory for Plotting
N = 1000;
t = linspace(t0, tf, N);

X_init = [X0; lambda0];
sol = ode45(@(t,X) EoM_minTime(t, X, B, mu_nd, umax), [t0 tf], X_init, opts);
X_eval = deval(sol, t);

pos = X_eval(1:3, :);       % nondimensional transfer position
vel = X_eval(4:6, :);
lambda = X_eval(7:12, :);

% dimensional transfer position for plotting
pos_dim = pos * lStar;

%% Calculate Control History
u = zeros(3, N);
for k = 1:N
    p = -B' * lambda(:,k);
    if norm(p) < 1e-10
        u(:,k) = zeros(3,1);
    else
        u(:,k) = (p / norm(p)) * umax;
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
yline(rEarth, 'Label','Earth Radius')
grid on;
xlabel('Time [hr]');
ylabel('Radius [km]');
title('Radius History');
legend('Transfer trajectory', 'Target evolution', 'Location', 'best');

%% Optional: Control magnitude vs Time
figure;
plot(t*tStar/3600, vecnorm(u,2,1), 'LineWidth', 1.5);
ylim([-1e-3, umax*1.1])
grid on;
xlabel('Time [hr]');
ylabel('||u|| [nondim]');
title('Control Magnitude History');

%% FUNCTION HELPER
function Xdot = twoBodyEOM_nd(~, X, mu_nd)
    r = X(1:3);
    v = X(4:6);

    r_norm = norm(r);

    a = -mu_nd * r / r_norm^3;

    Xdot = [v; a];
end