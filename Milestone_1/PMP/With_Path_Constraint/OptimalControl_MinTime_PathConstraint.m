clc; clear; close all;
% this script solves the minimum time optimal control transfer problem from
% parking orbit to target slot. It assumes no perturbations. Particularly,
% assumes that target orbit does not get purturbed.
%
% Earth intersection avoidance is enforced approximately through
% a running penalty in the cost functional, and rhoPenalty is increased
% gradually through homotopy continuation.
%
% IMPORTANT:
%   - kepler(M,e) returns E in degrees
%   - Parking_Orbit and Target store angles in degrees
%   - kep2cart is assumed to take angular inputs in degrees

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

opts = odeset('RelTol',1e-10, 'AbsTol',1e-10);

%% Earth avoidance / safety barrier settings
hSafe_km = 150;                    % minimum allowed altitude above Earth
hOn_km   = 300;                    % penalty turns on here
rSafe    = (rEarth + hSafe_km)/lStar;
rOn      = (rEarth + hOn_km)/lStar;
epsBar   = 1e-6;

%% Homotopy settings for rhoPenalty
rho_vec = [0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 5e-2];

%% Load Parking Orbit
load('Parking_Orbit.mat')

% get true anomaly in degrees:
E = kepler(Parking.M, Parking.ecc); 
ta = 2 * atan2(sqrt(1 + Parking.ecc) * sind(E/2), ...
               sqrt(1 - Parking.ecc) * cosd(E/2));
ta = rad2deg(ta);

% get initial cartesian state:
[x0, y0, z0, vx0, vy0, vz0] = kep2cart(Parking.a, Parking.ecc, Parking.inc, ...
    Parking.argp, Parking.raan, ta, mu);

% nondimensionalize:
r0 = [x0; y0; z0]/lStar;
v0 = [vx0; vy0; vz0]/vStar;
X0 = [r0;v0];

%% Solve Two-Point Boundary Value Problem (TPBVP)
% Optimization vector Z = [lambda1; lambda2; lambda3; lambda4; lambda5; lambda6; tf]
Z_guess = [4.7758, 5.6636, 3.3864, -2.9888, 5.5177, 7.4891, 5.9236]';
Z_guess = [1e-3; -1e-3; 1e-3; -1e-3; 1e-3; -1e-3; 10];
fsolve_opts = optimoptions('fsolve', ...
    'OptimalityTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'MaxFunctionEvaluations', 5000, ...
    'MaxIterations', 400, ...
    'Display', 'iter');

%% Homotopy loop on rhoPenalty
Z_sol = Z_guess;
exitflag_hist = nan(length(rho_vec),1);
fval_hist = nan(length(rho_vec),1);
rmin_hist = nan(length(rho_vec),1);
tf_hist = nan(length(rho_vec),1);

last_good_idx = 0;

for kRho = 1:length(rho_vec)

    rhoPenalty = rho_vec(kRho);

    fprintf('\n====================================================\n');
    fprintf('Homotopy step %d / %d : rhoPenalty = %.3e\n', ...
        kRho, length(rho_vec), rhoPenalty);
    fprintf('====================================================\n');

    [Z_try, fval, exitflag] = fsolve(@(Z) shooting_minTime_barrier(Z, X0, t0, ...
        umax, mu_nd, B, Target, opts, rSafe, rOn, rhoPenalty, epsBar), ...
        Z_sol, fsolve_opts);

    exitflag_hist(kRho) = exitflag;
    fval_hist(kRho) = norm(fval);

    if exitflag <= 0
        warning('Homotopy step failed at rhoPenalty = %.3e', rhoPenalty);
        break
    end

    % Accept solution and use as next initial guess
    Z_sol = Z_try;
    last_good_idx = kRho;

    % quick trajectory check
    lambda0_tmp = Z_sol(1:6);
    tf_tmp = Z_sol(7);
    tf_hist(kRho) = tf_tmp;

    X_init_tmp = [X0; lambda0_tmp];
    sol_tmp = ode45(@(t,X) EoM_minTime_barrier(t, X, B, mu_nd, umax, ...
        rSafe, rOn, rhoPenalty, epsBar), [t0 tf_tmp], X_init_tmp, opts);

    t_check = linspace(t0, tf_tmp, 400);
    X_check = deval(sol_tmp, t_check);
    r_check = X_check(1:3,:);
    r_norm_check = vecnorm(r_check,2,1);

    rmin_hist(kRho) = min(r_norm_check);

    fprintf('Converged at rhoPenalty = %.3e\n', rhoPenalty);
    fprintf('||shooting residual|| = %.3e\n', norm(fval));
    fprintf('tf = %.6f\n', tf_tmp);
    fprintf('minimum altitude = %.3f km\n', min(r_norm_check)*lStar - rEarth);
end

if last_good_idx == 0
    error('Homotopy failed at the first step. Try a better initial guess or smaller rhoPenalty values.')
end

%% Final converged solution
lambda0 = Z_sol(1:6);
tf = Z_sol(7);
rhoPenalty = rho_vec(last_good_idx);

fprintf('\nFinal accepted rhoPenalty: %.3e\n', rhoPenalty);
fprintf('Converged Initial Costate: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', lambda0);
fprintf('Converged Final Time (tf): %.4f\n', tf);

%% Generate Optimal Trajectory for Plotting
N = 1000;
t = linspace(t0, tf, N);

X_init = [X0; lambda0];
sol = ode45(@(t,X) EoM_minTime_barrier(t, X, B, mu_nd, umax, ...
    rSafe, rOn, rhoPenalty, epsBar), [t0 tf], X_init, opts);
X_eval = deval(sol, t);

pos = X_eval(1:3, :);       % nondimensional transfer position
vel = X_eval(4:6, :);
lambda = X_eval(7:12, :);

% dimensional transfer position for plotting
pos_dim = pos * lStar;

%% Calculate Control History
u = zeros(3, N);
phi_hist = zeros(1, N);
rmag_hist = zeros(1, N);
alt_hist = zeros(1, N);

for k = 1:N
    p = -B' * lambda(:,k);

    if norm(p) < 1e-10
        u(:,k) = zeros(3,1);
    else
        u(:,k) = (p / norm(p)) * umax;
    end

    rmag_hist(k) = norm(pos(:,k)) * lStar;
    alt_hist(k) = rmag_hist(k) - rEarth;

    [phi_hist(k), ~] = radiusPenaltyCost(pos(:,k), rSafe, rOn, epsBar);
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
    % Target.M and Target.n_nd usage remains consistent:
    % Target.M is in degrees, Target.n_nd is rad / nondim-time
    M_k = Target.M + rad2deg(Target.n_nd * (t(k) - t0));
    M_k = mod(M_k, 360);

    % kepler returns E in degrees
    E_k = kepler(M_k, Target.ecc);

    % convert to true anomaly in degrees
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

% safe radius sphere
surf((rSafe*lStar)*xe, (rSafe*lStar)*ye, (rSafe*lStar)*ze, ...
    'FaceColor', [1.0 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.08);

plot3(rPark_dim(1,1), rPark_dim(2,1), rPark_dim(3,1), 'ko', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 6);
plot3(rTarget_dim(1,1), rTarget_dim(2,1), rTarget_dim(3,1), 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 6);
plot3(pos_dim(1,end), pos_dim(2,end), pos_dim(3,end), 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 6);

legend('Earth', 'Initial parking orbit', 'Target state evolution', ...
       'Optimal transfer arc', 'Safe radius', ...
       'Initial spacecraft state', 'Initial target state', ...
       'Final transfer state', 'Location', 'best');

view(3);

%% Optional: Position vs Time
figure;
plot(t*tStar/3600, vecnorm(pos_dim,2,1), 'b', 'LineWidth', 1.5); hold on;
plot(t*tStar/3600, vecnorm(rTarget_dim,2,1), 'r--', 'LineWidth', 1.5);
yline(rSafe*lStar, 'k--', 'LineWidth', 1.2);
yline(rOn*lStar, 'm--', 'LineWidth', 1.2);
grid on;
xlabel('Time [hr]');
ylabel('Radius [km]');
title('Radius History');
legend('Transfer trajectory', 'Target evolution', 'r_{safe}', 'r_{on}', 'Location', 'best');

%% Optional: Altitude vs Time
figure;
plot(t*tStar/3600, alt_hist, 'b', 'LineWidth', 1.5); hold on;
yline(hSafe_km, 'k--', 'LineWidth', 1.2);
yline(hOn_km, 'm--', 'LineWidth', 1.2);
grid on;
xlabel('Time [hr]');
ylabel('Altitude [km]');
title('Altitude History');
legend('Transfer trajectory', 'h_{safe}', 'h_{on}', 'Location', 'best');

%% Optional: Control magnitude vs Time
figure;
plot(t*tStar/3600, vecnorm(u,2,1), 'LineWidth', 1.5);
grid on;
xlabel('Time [hr]');
ylabel('||u|| [nondim]');
title('Control Magnitude History');

%% Optional: Penalty vs Time
figure;
plot(t*tStar/3600, phi_hist, 'LineWidth', 1.5);
grid on;
xlabel('Time [hr]');
ylabel('\Phi(r)');
title('Earth-Avoidance Running Penalty');

%% Optional: Homotopy diagnostics
valid_idx = ~isnan(exitflag_hist);

figure;
subplot(2,1,1)
semilogx(rho_vec(valid_idx), fval_hist(valid_idx), 'o-', 'LineWidth', 1.5)
grid on
xlabel('\rho_{penalty}')
ylabel('||\Psi||')
title('Homotopy Residual History')

subplot(2,1,2)
semilogx(rho_vec(valid_idx), rmin_hist(valid_idx)*lStar - rEarth, 'o-', 'LineWidth', 1.5)
grid on
xlabel('\rho_{penalty}')
ylabel('Minimum Altitude [km]')
title('Homotopy Minimum Altitude History')

%% FUNCTION HELPER
function Xdot = twoBodyEOM_nd(~, X, mu_nd)
    r = X(1:3);
    v = X(4:6);

    r_norm = norm(r);

    a = -mu_nd * r / r_norm^3;

    Xdot = [v; a];
end