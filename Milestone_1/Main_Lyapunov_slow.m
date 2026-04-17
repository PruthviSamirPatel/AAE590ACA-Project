clc; clear; close all;

fprintf('\n============================================================\n');
fprintf(' Simple Slow-Variable Lyapunov Orbit Transfer\n');
fprintf('============================================================\n');

%% Load environment
load('Earth_params.mat');
rEarth = Earth.radius;      % km
mu_km  = Earth.mu;          % km^3/s^2

%% Load initial and target orbits
load('Parking_Orbit.mat');
load('Orbital_Slot.mat');

%% Canonical scaling
LU = Parking.a;                 % use initial semimajor axis as length unit
TU = sqrt(LU^3 / mu_km);
VU = LU / TU;
AU = LU / TU^2;
mu = 1.0;                       % nondimensional

Re_nd = rEarth / LU;

fprintf('\nScaling:\n');
fprintf('  LU = %.6f km\n', LU);
fprintf('  TU = %.6f s\n', TU);
fprintf('  VU = %.6f km/s\n', VU);
fprintf('  AU = %.9e km/s^2\n', AU);

%% Control authority
umax_dim = 1e-5;                % km/s^2
umax_nd  = umax_dim / AU;

fprintf('\nControl:\n');
fprintf('  umax_dim = %.6e km/s^2\n', umax_dim);
fprintf('  umax_nd  = %.6e\n', umax_nd);

%% Initial / target Keplerian states (nondim a, radians otherwise)
x0 = [ ...
    Parking.a / LU;
    Parking.ecc;
    deg2rad(Parking.inc);
    deg2rad(Parking.raan);
    deg2rad(Parking.argp);
    deg2rad(Parking.M)];

xT = [ ...
    Target.a / LU;
    Target.ecc;
    deg2rad(Target.inc);
    deg2rad(Target.raan);
    deg2rad(Target.argp);
    deg2rad(Target.M)];

fprintf('\nInitial orbit:\n');
fprintf('  a0 = %.6f km, e0 = %.6f, i0 = %.6f deg, RAAN0 = %.6f deg, omega0 = %.6f deg, M0 = %.6f deg\n', ...
    Parking.a, Parking.ecc, Parking.inc, Parking.raan, Parking.argp, Parking.M);

fprintf('\nTarget orbit:\n');
fprintf('  aT = %.6f km, eT = %.6f, iT = %.6f deg, RAANT = %.6f deg, omegaT = %.6f deg, MT = %.6f deg\n', ...
    Target.a, Target.ecc, Target.inc, Target.raan, Target.argp, Target.M);

%% Lyapunov matrices
K = eye(5);
P = eye(5);

%% Integration settings
tf_days = 7.0;
tf_sec  = tf_days * 24 * 3600;
tf_nd   = tf_sec / TU;

dt_nd = 0.002;                        % fixed RK4 step
N = floor(tf_nd / dt_nd) + 1;
t_nd = linspace(0, dt_nd*(N-1), N).';
t_days = t_nd * TU / 86400;

fprintf('\nIntegration settings:\n');
fprintf('  tf = %.6f days\n', tf_days);
fprintf('  dt = %.6e nondim\n', dt_nd);
fprintf('  N  = %d\n', N);

%% Penalty settings
Pen.Re_nd   = Re_nd;
Pen.r_warn  = (rEarth + 300)/LU;
Pen.r_min   = (rEarth + 150)/LU;
Pen.k_pen   = 0.02;

%% Simulate with saturation
xHist       = zeros(N,6);
uHist       = zeros(N,3);
uMagHist    = zeros(N,1);
dxslowHist  = zeros(N,5);
rEciHist_nd = zeros(N,3);
penHist     = zeros(N,1);

xHist(1,:) = x0.';

for k = 1:N-1
    xk = xHist(k,:).';

    [k1, aux1] = dyn_kepler_controller_sat_penalty(xk, xT, mu, K, P, true, umax_nd, Pen);
    [k2, ~]    = dyn_kepler_controller_sat_penalty(xk + 0.5*dt_nd*k1, xT, mu, K, P, true, umax_nd, Pen);
    [k3, ~]    = dyn_kepler_controller_sat_penalty(xk + 0.5*dt_nd*k2, xT, mu, K, P, true, umax_nd, Pen);
    [k4, ~]    = dyn_kepler_controller_sat_penalty(xk + dt_nd*k3,     xT, mu, K, P, true, umax_nd, Pen);

    xnext = xk + (dt_nd/6)*(k1 + 2*k2 + 2*k3 + k4);

    % guards
    xnext(1) = max(xnext(1), 1e-4);
    xnext(2) = min(max(xnext(2), 1e-6), 0.95);
    xnext(3) = min(max(xnext(3), 1e-4), pi-1e-4);
    xnext(4) = mod(xnext(4), 2*pi);
    xnext(5) = mod(xnext(5), 2*pi);
    xnext(6) = mod(xnext(6), 2*pi);

    xHist(k+1,:) = xnext.';

    uHist(k,:)       = aux1.u.';
    uMagHist(k)      = norm(aux1.u);
    dxslowHist(k,:)  = aux1.dxslow.';
    rEciHist_nd(k,:) = aux1.r_eci.';
    penHist(k)       = aux1.penalty;
end

[~, auxN] = dyn_kepler_controller_sat_penalty(xHist(end,:).', xT, mu, K, P, true, umax_nd, Pen);
uHist(end,:)       = auxN.u.';
uMagHist(end)      = norm(auxN.u);
dxslowHist(end,:)  = auxN.dxslow.';
rEciHist_nd(end,:) = auxN.r_eci.';
penHist(end)       = auxN.penalty;

%% Convert to physical units
aHist_km      = xHist(:,1) * LU;
eHist         = xHist(:,2);
iHist_deg     = rad2deg(xHist(:,3));
OmegaHist_deg = rad2deg(xHist(:,4));
omegaHist_deg = rad2deg(xHist(:,5));
MHist_deg     = rad2deg(xHist(:,6));

uMag_kms2 = uMagHist * AU;
rEciHist_km = rEciHist_nd * LU;

%% Safety metrics
altHist_km = vecnorm(rEciHist_km,2,2) - rEarth;
minAlt_km  = min(altHist_km);

%% Final errors
dx_final = dxslowHist(end,:);

fprintf('\n============================================================\n');
fprintf(' RESULTS\n');
fprintf('============================================================\n');
fprintf('Final delta a/LU    = %+12.6e\n', dx_final(1));
fprintf('Final delta e       = %+12.6e\n', dx_final(2));
fprintf('Final delta i       = %+12.6e rad\n', dx_final(3));
fprintf('Final delta Omega   = %+12.6e rad\n', dx_final(4));
fprintf('Final delta omega   = %+12.6e rad\n', dx_final(5));
fprintf('Peak ||u||          = %.6e km/s^2 = %.6f m/s^2\n', ...
        max(uMag_kms2), max(uMag_kms2)*1e3);
fprintf('Approx delta-v      = %.6f km/s\n', trapz(t_nd, uMagHist) * VU);
fprintf('Minimum altitude    = %.6f km\n', minAlt_km);

%% 3D plot
figure('Color','w','Position',[100 100 900 650]);
hold on; grid on; axis equal;

% Earth
[xe, ye, ze] = sphere(80);
surf(rEarth*xe, rEarth*ye, rEarth*ze, ...
    'FaceAlpha',0.15,'EdgeColor','none');

% Transfer trajectory
plot3(rEciHist_km(:,1), rEciHist_km(:,2), rEciHist_km(:,3), ...
    'b', 'LineWidth', 2);

% Start and end markers
plot3(rEciHist_km(1,1), rEciHist_km(1,2), rEciHist_km(1,3), ...
    'go','MarkerFaceColor','g','MarkerSize',8);

plot3(rEciHist_km(end,1), rEciHist_km(end,2), rEciHist_km(end,3), ...
    'ro','MarkerFaceColor','r','MarkerSize',8);

% Axis labels
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');

title('3D Orbit Transfer (Slow-Variable Lyapunov)');

legend('Earth','Transfer Trajectory','Start','End','Location','best');

% Force 3D view
view(35,25);

% Improve lighting
camlight; lighting gouraud;

figure;
plot(t_days, uMag_kms2*1e3, 'LineWidth',1.5); grid on;
xlabel('t [days]');
ylabel('||u|| [m/s^2]');
title('Control Magnitude');

figure;
subplot(5,1,1); plot(t_days, dxslowHist(:,1), 'LineWidth',1.2); grid on; ylabel('\delta a/LU');
subplot(5,1,2); plot(t_days, dxslowHist(:,2), 'LineWidth',1.2); grid on; ylabel('\delta e');
subplot(5,1,3); plot(t_days, rad2deg(dxslowHist(:,3)), 'LineWidth',1.2); grid on; ylabel('\delta i [deg]');
subplot(5,1,4); plot(t_days, rad2deg(dxslowHist(:,4)), 'LineWidth',1.2); grid on; ylabel('\delta \Omega [deg]');
subplot(5,1,5); plot(t_days, rad2deg(dxslowHist(:,5)), 'LineWidth',1.2); grid on; ylabel('\delta \omega [deg]');
xlabel('t [days]');
sgtitle('Slow-Variable Errors');

figure;
plot(t_days, altHist_km, 'LineWidth',1.5); hold on;
yline(150,'--r','h_{min}');
yline(300,'--m','h_{warn}');
grid on;
xlabel('t [days]');
ylabel('Altitude [km]');
title('Altitude History');

%% Save
RefSlow.t_days      = t_days;
RefSlow.xHist       = xHist;
RefSlow.dxslowHist  = dxslowHist;
RefSlow.uHist       = uHist;
RefSlow.uMag_kms2   = uMag_kms2;
RefSlow.rEciHist_km = rEciHist_km;
RefSlow.altHist_km  = altHist_km;
RefSlow.penHist     = penHist;
save('Reference_Trajectory_Lyapunov_Slow.mat','RefSlow');

fprintf('\nSaved: Reference_Trajectory_Lyapunov_Slow.mat\n');
fprintf('Done.\n');

%% ========================================================================
function [xdot, aux] = dyn_kepler_controller_sat_penalty(x, xT, mu, K, P, useSat, umax, Pen)
    % x = [a; e; i; Omega; omega; M]

    a   = max(x(1),1e-4);
    e   = min(max(x(2),1e-6),0.95);
    inc = min(max(x(3),1e-4),pi-1e-4);
    Om  = x(4);
    w   = x(5);
    M   = x(6);

    E = kepler_fast(M, e);

    sinf = sqrt(1-e^2) * sin(E) / (1 - e*cos(E));
    cosf = (cos(E) - e) / (1 - e*cos(E));
    f = atan2(sinf, cosf);

    p = a*(1-e^2);
    h = sqrt(mu*p);
    r = p/(1 + e*cos(f));

    % Slow-variable Gauss matrix
    B = zeros(5,3);

    B(1,1) = 2*a^2/h * e*sin(f);
    B(1,2) = 2*a^2/h * (p/r);

    B(2,1) = p*sin(f)/h;
    B(2,2) = ((p+r)*cos(f) + r*e)/h;

    B(3,3) = r*cos(w+f)/h;
    B(4,3) = r*sin(w+f)/(h*sin(inc));

    % guard small e
    ee = max(e,1e-4);
    B(5,1) = -(p*cos(f))/(h*ee);
    B(5,2) =  (p+r)*sin(f)/(h*ee);
    B(5,3) = -(r*sin(w+f)*cos(inc))/(h*sin(inc));

    dxslow = x(1:5) - xT(1:5);

    H = B.' * K.' * K * B;
    g = B.' * K.' * P * dxslow;

    u_nom = -(H + 1e-8*eye(3)) \ g;

    % simple radial penalty in RTN coordinates
    if r < Pen.r_warn
        sigma = (Pen.r_warn - r) / (Pen.r_warn - Pen.r_min);
        sigma = max(0, sigma);
        u_pen = [Pen.k_pen * sigma^2; 0; 0];
    else
        u_pen = [0;0;0];
    end

    u_cmd = u_nom + u_pen;

    if useSat
        nu = norm(u_cmd);
        if nu > umax
            u = (umax/nu) * u_cmd;
        else
            u = u_cmd;
        end
    else
        u = u_cmd;
    end

    n = sqrt(mu / a^3);

    xdot = zeros(6,1);
    xdot(1:5) = B*u;
    xdot(6)   = n;

    % ECI position for plotting
    r_pf = [r*cos(f); r*sin(f); 0];

    R3_Om = [ cos(Om) -sin(Om) 0;
              sin(Om)  cos(Om) 0;
              0        0       1 ];
    R1_i  = [ 1 0 0;
              0 cos(inc) -sin(inc);
              0 sin(inc)  cos(inc) ];
    R3_w  = [ cos(w) -sin(w) 0;
              sin(w)  cos(w) 0;
              0       0      1 ];

    r_eci = R3_Om * R1_i * R3_w * r_pf;

    aux.u = u;
    aux.dxslow = dxslow;
    aux.r_eci = r_eci;
    aux.penalty = norm(u_pen);
end

function E = kepler_fast(M, e)
    M = mod(M + pi, 2*pi) - pi;

    if e < 0.8
        E = M;
    else
        E = pi;
    end

    for k = 1:5
        F  = E - e*sin(E) - M;
        dF = 1 - e*cos(E);
        E  = E - F/dF;
    end
end