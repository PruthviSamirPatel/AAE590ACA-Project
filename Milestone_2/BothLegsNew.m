clc; clear; close all

%% ============================================================
% SIMPLE TWO-PHASE LOW-THRUST TRANSFER
%
% Phase 1:
%   - transfer from parking orbit to a drift orbit
%   - drift orbit is chosen so J2 RAAN drift matches target after coasting
%   - then coast analytically until RAAN match
%
% Phase 2:
%   - Lyapunov transfer from drift orbit to target
%   - uses a fixed target perigee during the transfer
%
% Features:
%   - fixed-step RK4 (predictable runtime)
%   - J2 included in active propagation
%   - Earth-radius / safe-perigee constraints
%   - simple MEE Lyapunov control
%
% Units:
%   km, s, rad
%% ============================================================

%% ------------------------- EARTH ----------------------------
Earth.mu     = 398600.4418;         % km^3/s^2
Earth.radius = 6378.1363;           % km
Earth.J2     = 1.08262668e-3;

mu     = Earth.mu;
rEarth = Earth.radius;
J2     = Earth.J2;

%% ---------------------- USER SETTINGS -----------------------
% Parking orbit
Parking.alt  = 610;                 % km
Parking.a    = rEarth + Parking.alt;
Parking.ecc  = 0.001;
Parking.inc  = deg2rad(51.9);
Parking.raan = deg2rad(60.0);
Parking.argp = deg2rad(0.0);
Parking.M    = deg2rad(10.0);

% Target orbit
Target.alt   = 630;                 % km
Target.a     = rEarth + Target.alt;
Target.ecc   = 0.0005;
Target.inc   = deg2rad(51.9);
Target.raan  = deg2rad(40.0);
Target.argp  = deg2rad(0.0);
Target.M     = deg2rad(30.0);

% Transfer settings
cfg.dt                = 30.0;       % fixed RK4 step [s]
cfg.uMax_km_s2        = 1e-5;       % max accel [km/s^2] = 0.01 m/s^2
cfg.hSafe_km          = 150;        % hard minimum altitude
cfg.hGuard_km         = 220;        % start safety recovery below this perigee
cfg.hBarrier_km       = 300;        % current-radius outward push below this

cfg.phase1_max_days   = 1.5;        % max active time for transfer to drift orbit
cfg.coast_nom_days    = 40.0;       % nominal coast time used to size drift orbit
cfg.phase2_max_days   = 4.0;        % max active time for final transfer

cfg.tol_p_km          = 5.0;        % phase 1 convergence on p
cfg.tol_e             = 2e-3;       % phase 1 convergence on e
cfg.rTol_km           = 1.0;        % phase 2 rendezvous tolerance
cfg.vTol_km_s         = 1e-3;

% Lyapunov gains
cfg.k1                = 2e-2;       % phase 1 gain
cfg.k2                = 3e-2;       % phase 2 gain
cfg.kSafe             = 4e-2;       % safety recovery gain
cfg.reg               = 1e-8;

%% ---------------- INITIAL STATES ----------------------------
nuC0 = mean2true(Parking.M, Parking.ecc);
nuT0 = mean2true(Target.M,  Target.ecc);

[rC0, vC0] = kep2rv(Parking.a, Parking.ecc, Parking.inc, ...
                    Parking.raan, Parking.argp, nuC0, mu);

[rT0, vT0] = kep2rv(Target.a, Target.ecc, Target.inc, ...
                    Target.raan, Target.argp, nuT0, mu);

xC = [rC0; vC0];
xT = [rT0; vT0];

%% ---------------- DRIFT ORBIT SELECTION ---------------------
deltaOm0 = angwrap(Parking.raan - Target.raan);  % chaser - target
raanDotT = secular_raan_rate(Target.a, Target.ecc, Target.inc, Earth);

% Need: deltaOm0 + (raanDotD - raanDotT)*t = 0
raanDotReq = raanDotT - deltaOm0/(cfg.coast_nom_days*86400);

eDrift = 1e-3;
iDrift = Parking.inc;

aMin = rEarth + max(cfg.hSafe_km + 50, 200);
aMax = Target.a - 1.0;

f = @(a) secular_raan_rate(a, eDrift, iDrift, Earth) - raanDotReq;

if f(aMin)*f(aMax) > 0
    error('Could not bracket a feasible drift orbit. Increase coast time or lower safety floor.');
end

aDrift = fzero(f, [aMin, aMax]);
pDrift = aDrift*(1 - eDrift^2);

fprintf('\n================ DRIFT ORBIT ========================\n');
fprintf('Nominal coast time       = %.3f days\n', cfg.coast_nom_days);
fprintf('Required drift SMA       = %.3f km\n', aDrift);
fprintf('Required drift altitude  = %.3f km\n', aDrift - rEarth);
fprintf('=====================================================\n\n');

%% ---------------- LOGGING SETUP -----------------------------
Tlog  = [];
XClog = [];
XTlog = [];
Ulog  = [];
Phase = [];

tNow = 0.0;

%% ============================================================
% PHASE 1: ACTIVE TRANSFER TO DRIFT ORBIT
%% ============================================================
N1 = ceil(cfg.phase1_max_days*86400/cfg.dt);
done1 = false;

for k = 1:N1
    meeC = rv2mee(xC, mu);

    pC = meeC(1);
    eC = hypot(meeC(2), meeC(3));

    % Desired drift orbit: same plane, near-circular
    meeDes1 = [pDrift; ...
               eDrift*cos(Parking.raan + Parking.argp); ...
               eDrift*sin(Parking.raan + Parking.argp); ...
               meeC(4); ...
               meeC(5); ...
               meeC(6)];

    uC = control_phase1(xC, meeDes1, cfg, Earth);

    % Propagate chaser and target
    xC = rk4_step(@(x) rhs_orbit(x, uC, Earth), xC, cfg.dt);
    xT = rk4_step(@(x) rhs_orbit(x, [0;0;0], Earth), xT, cfg.dt);

    tNow = tNow + cfg.dt;

    % Log
    Tlog  = [Tlog; tNow]; %#ok<AGROW>
    XClog = [XClog; xC.']; %#ok<AGROW>
    XTlog = [XTlog; xT.']; %#ok<AGROW>
    Ulog  = [Ulog;  uC.']; %#ok<AGROW>
    Phase = [Phase; 1]; %#ok<AGROW>

    % Safety stop
    [rpAlt, rAlt] = altitude_metrics(xC, Earth);
    if rpAlt < cfg.hSafe_km || rAlt < cfg.hSafe_km
        error('Phase 1 violated Earth safety constraint.');
    end

    % Convergence to drift orbit
    if abs(pC - pDrift) < cfg.tol_p_km && eC < cfg.tol_e
        done1 = true;
        fprintf('Phase 1 converged in %.3f hours.\n', tNow/3600);
        break
    end
end

if ~done1
    fprintf('Phase 1 hit max time without full convergence. Continuing anyway.\n');
end

%% ============================================================
% PHASE 1 COAST: ANALYTIC J2 DRIFT UNTIL RAAN MATCH
%% ============================================================
coeC1 = rv2coe(xC, mu);
coeT1 = rv2coe(xT, mu);

raanDotC = secular_raan_rate(coeC1.a, coeC1.e, coeC1.i, Earth);
raanDotT = secular_raan_rate(coeT1.a, coeT1.e, coeT1.i, Earth);

deltaOm1 = angwrap(coeC1.raan - coeT1.raan);
dRate    = raanDotC - raanDotT;

if abs(dRate) < 1e-14
    error('Relative RAAN drift is too small. Change drift orbit or coast time.');
end

tCoast = -deltaOm1 / dRate;

if tCoast < 0
    error('Computed coast time is negative. Chosen drift orbit does not close the RAAN error.');
end

fprintf('Exact coast time to RAAN match = %.6f days\n', tCoast/86400);

coeC2 = coast_orbit_J2(coeC1, tCoast, Earth);
coeT2 = coast_orbit_J2(coeT1, tCoast, Earth);

[xC, xT] = deal([kep2rv_state(coeC2, mu)], [kep2rv_state(coeT2, mu)]);

tNow = tNow + tCoast;

% log coast endpoints only
Tlog  = [Tlog; tNow]; %#ok<AGROW>
XClog = [XClog; xC.']; %#ok<AGROW>
XTlog = [XTlog; xT.']; %#ok<AGROW>
Ulog  = [Ulog; [0 0 0]]; %#ok<AGROW>
Phase = [Phase; 1]; %#ok<AGROW>

fprintf('RAAN error after coast = %.6e deg\n', rad2deg(angwrap(coeC2.raan - coeT2.raan)));

%% ============================================================
% PHASE 2: FIXED-PERIGEE LYAPUNOV TRANSFER TO TARGET
%% ============================================================
coeC2s = rv2coe(xC, mu);
rpFix  = Target.a*(1 - Target.ecc);        % fixed target perigee
ra0    = coeC2s.a*(1 + coeC2s.e);          % start apogee
raT    = Target.a*(1 + Target.ecc);        % target apogee

N2 = ceil(cfg.phase2_max_days*86400/cfg.dt);
done2 = false;

for k = 1:N2
    tau = (k-1)/(N2-1);
    raDes = (1-tau)*ra0 + tau*raT;
    aDes  = 0.5*(rpFix + raDes);
    eDes  = (raDes - rpFix)/(raDes + rpFix);
    pDes  = aDes*(1 - eDes^2);

    coeT = rv2coe(xT, mu);

    lonPerDes = coeT.raan + coeT.argp;
    fDes = eDes*cos(lonPerDes);
    gDes = eDes*sin(lonPerDes);
    hDes = tan(coeT.i/2)*cos(coeT.raan);
    kDes = tan(coeT.i/2)*sin(coeT.raan);
    LDes = wrap2pi(coeT.raan + coeT.argp + coeT.nu);

    meeDes2 = [pDes; fDes; gDes; hDes; kDes; LDes];

    uC = control_phase2(xC, meeDes2, cfg, Earth);

    xC = rk4_step(@(x) rhs_orbit(x, uC, Earth), xC, cfg.dt);
    xT = rk4_step(@(x) rhs_orbit(x, [0;0;0], Earth), xT, cfg.dt);

    tNow = tNow + cfg.dt;

    Tlog  = [Tlog; tNow]; %#ok<AGROW>
    XClog = [XClog; xC.']; %#ok<AGROW>
    XTlog = [XTlog; xT.']; %#ok<AGROW>
    Ulog  = [Ulog;  uC.']; %#ok<AGROW>
    Phase = [Phase; 2]; %#ok<AGROW>

    [rpAlt, rAlt] = altitude_metrics(xC, Earth);
    if rpAlt < cfg.hSafe_km || rAlt < cfg.hSafe_km
        error('Phase 2 violated Earth safety constraint.');
    end

    dr = norm(xC(1:3) - xT(1:3));
    dv = norm(xC(4:6) - xT(4:6));

    if dr < cfg.rTol_km && dv < cfg.vTol_km_s
        done2 = true;
        fprintf('Phase 2 rendezvous reached.\n');
        break
    end
end

if ~done2
    fprintf('Phase 2 ended at max time.\n');
end

%% --------------------- POST PROCESS -------------------------
N = size(XClog,1);
tDays = Tlog/86400;

rC = XClog(:,1:3);
vC = XClog(:,4:6);
rT = XTlog(:,1:3);
vT = XTlog(:,4:6);

aC = zeros(N,1); eC = zeros(N,1); iC = zeros(N,1); OmC = zeros(N,1);
aT = zeros(N,1); eT = zeros(N,1); iT = zeros(N,1); OmT = zeros(N,1);

rpAlt = zeros(N,1);
rAlt  = zeros(N,1);
drHist = zeros(N,1);
dvHist = zeros(N,1);
raanErr = zeros(N,1);

for k = 1:N
    coeC = rv2coe(XClog(k,:).', mu);
    coeT = rv2coe(XTlog(k,:).', mu);

    aC(k) = coeC.a; eC(k) = coeC.e; iC(k) = rad2deg(coeC.i); OmC(k) = rad2deg(coeC.raan);
    aT(k) = coeT.a; eT(k) = coeT.e; iT(k) = rad2deg(coeT.i); OmT(k) = rad2deg(coeT.raan);

    [rpAlt(k), rAlt(k)] = altitude_metrics(XClog(k,:).', Earth);
    drHist(k) = norm(rC(k,:) - rT(k,:));
    dvHist(k) = norm(vC(k,:) - vT(k,:));
    raanErr(k) = rad2deg(angwrap(coeC.raan - coeT.raan));
end

uNorm = vecnorm(Ulog,2,2);
DeltaV = trapz(Tlog, uNorm);

fprintf('\n================ FINAL SUMMARY ======================\n');
fprintf('Total elapsed time         = %.6f days\n', tDays(end));
fprintf('Integrated control DV      = %.6f km/s\n', DeltaV);
fprintf('Final ||dr||               = %.6e km\n', drHist(end));
fprintf('Final ||dv||               = %.6e km/s\n', dvHist(end));
fprintf('Final RAAN error           = %.6e deg\n', raanErr(end));
fprintf('Final perigee altitude     = %.6f km\n', rpAlt(end));
fprintf('=====================================================\n\n');

%% ------------------------- PLOTS ----------------------------
figure('Color','w');
plot3(rC(:,1), rC(:,2), rC(:,3), 'r', 'LineWidth',1.2); hold on
plot3(rT(:,1), rT(:,2), rT(:,3), 'b--', 'LineWidth',1.2);
[sx,sy,sz] = sphere(50);
surf(rEarth*sx, rEarth*sy, rEarth*sz, 'FaceAlpha',0.15, 'EdgeColor','none');
grid on; axis equal; view(3)
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Two-Phase Transfer');

figure('Color','w');
plot(tDays, raanErr, 'LineWidth',1.5); grid on
xlabel('t [days]'); ylabel('\Omega_C - \Omega_T [deg]');
title('RAAN Error');

figure('Color','w');
plot(tDays, aC, 'LineWidth',1.5); hold on
plot(tDays, aT, '--', 'LineWidth',1.5);
grid on
xlabel('t [days]'); ylabel('a [km]');
title('Semimajor Axis');

figure('Color','w');
plot(tDays, rpAlt, 'LineWidth',1.5); hold on
plot(tDays, rAlt, 'LineWidth',1.2);
yline(cfg.hSafe_km, 'r--', 'LineWidth',1.3);
yline(cfg.hGuard_km, 'm--', 'LineWidth',1.0);
grid on
xlabel('t [days]'); ylabel('Altitude [km]');
title('Perigee / Current Altitude');
legend('Perigee altitude','Current altitude','h_{safe}','h_{guard}','Location','best');

figure('Color','w');
plot(tDays, drHist, 'LineWidth',1.5); hold on
plot(tDays, dvHist, 'LineWidth',1.5);
grid on
xlabel('t [days]'); ylabel('Norm');
title('Relative State');
legend('||\Delta r|| [km]','||\Delta v|| [km/s]','Location','best');

figure('Color','w');
plot(tDays, 1000*Ulog(:,1), 'LineWidth',1.0); hold on
plot(tDays, 1000*Ulog(:,2), 'LineWidth',1.0);
plot(tDays, 1000*Ulog(:,3), 'LineWidth',1.0);
plot(tDays, 1000*uNorm, 'k', 'LineWidth',1.5);
grid on
xlabel('t [days]'); ylabel('u [m/s^2]');
title(sprintf('RTN Control History, \\Delta v = %.4f km/s', DeltaV));
legend('u_R','u_T','u_N','||u||','Location','best');

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function u = control_phase1(x, meeDes, cfg, Earth)
    mu = Earth.mu;

    mee = rv2mee(x, mu);
    B   = mee_B(mee, mu);

    p = mee(1); f = mee(2); g = mee(3);

    pDes = meeDes(1); fDes = meeDes(2); gDes = meeDes(3);

    evec = [(p - pDes)/max(pDes,1);
            f - fDes;
            g - gDes];

    A = B(1:3,1:2);
    W = diag([1, 10, 10]);

    u12 = -cfg.k1 * ((A'*W*A + cfg.reg*eye(2)) \ (A'*W*evec));
    u   = [u12; 0];

    u = add_safety_override(x, u, cfg, Earth);

    nu = norm(u);
    if nu > cfg.uMax_km_s2
        u = cfg.uMax_km_s2 * u / nu;
    end
end

function u = control_phase2(x, meeDes, cfg, Earth)
    mu = Earth.mu;

    mee = rv2mee(x, mu);
    B   = mee_B(mee, mu);

    evec = [ (mee(1) - meeDes(1))/max(meeDes(1),1);
             mee(2) - meeDes(2);
             mee(3) - meeDes(3);
             mee(4) - meeDes(4);
             mee(5) - meeDes(5);
             angwrap(mee(6) - meeDes(6)) ];

    W = diag([1, 8, 8, 4, 4, 0.2]);

    u = -cfg.k2 * ((B'*W*B + cfg.reg*eye(3)) \ (B'*W*evec));

    u = add_safety_override(x, u, cfg, Earth);

    nu = norm(u);
    if nu > cfg.uMax_km_s2
        u = cfg.uMax_km_s2 * u / nu;
    end
end

function u = add_safety_override(x, uIn, cfg, Earth)
    mu = Earth.mu;

    [rpAlt, rAlt] = altitude_metrics(x, Earth);

    if rpAlt >= cfg.hGuard_km
        u = uIn;

        % light outward push if current altitude gets low
        if rAlt < cfg.hBarrier_km
            u(1) = u(1) + 2e-6;
        end
        return
    end

    % Safety recovery: raise perigee and damp eccentricity
    mee = rv2mee(x, mu);
    B   = mee_B(mee, mu);

    p = mee(1);
    f = mee(2);
    g = mee(3);

    pSafe = (Earth.radius + cfg.hGuard_km) * (1 + hypot(f,g));

    evec = [ (p - pSafe)/max(pSafe,1);
             f;
             g ];

    A = B(1:3,1:2);
    W = diag([1, 10, 10]);

    u12 = -cfg.kSafe * ((A'*W*A + cfg.reg*eye(2)) \ (A'*W*evec));
    u   = [u12; 0];

    nu = norm(u);
    if nu > cfg.uMax_km_s2
        u = cfg.uMax_km_s2 * u / nu;
    end
end

function xnext = rk4_step(fun, x, dt)
    k1 = fun(x);
    k2 = fun(x + 0.5*dt*k1);
    k3 = fun(x + 0.5*dt*k2);
    k4 = fun(x + dt*k3);
    xnext = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
end

function xdot = rhs_orbit(x, uRTN, Earth)
    r = x(1:3);
    v = x(4:6);

    mu = Earth.mu;
    J2 = Earth.J2;
    Re = Earth.radius;

    rmag = norm(r);

    a2b = -mu * r / rmag^3;
    aJ2 = accelJ2(r, mu, J2, Re);
    aU  = rtn_to_eci(r, v, uRTN);

    xdot = [v; a2b + aJ2 + aU];
end

function aJ2 = accelJ2(r, mu, J2, Re)
    x = r(1); y = r(2); z = r(3);
    r2 = dot(r,r);
    r1 = sqrt(r2);
    z2 = z^2;

    fac = 1.5 * J2 * mu * Re^2 / r1^5;
    cxy = 5*z2/r2 - 1;
    cz  = 5*z2/r2 - 3;

    aJ2 = fac * [x*cxy; y*cxy; z*cz];
end

function aECI = rtn_to_eci(r, v, aRTN)
    Rhat = r / norm(r);
    hvec = cross(r, v);
    Nhat = hvec / norm(hvec);
    That = cross(Nhat, Rhat);

    Q = [Rhat, That, Nhat];
    aECI = Q * aRTN;
end

function [rpAlt, rAlt] = altitude_metrics(x, Earth)
    coe = rv2coe(x, Earth.mu);
    rp  = coe.a*(1 - coe.e);

    rpAlt = rp - Earth.radius;
    rAlt  = norm(x(1:3)) - Earth.radius;
end

function rate = secular_raan_rate(a, e, inc, Earth)
    mu = Earth.mu; J2 = Earth.J2; Re = Earth.radius;
    p = a*(1 - e^2);
    n = sqrt(mu/a^3);
    rate = -(3/2)*J2*n*(Re/p)^2*cos(inc);
end

function rate = secular_argp_rate(a, e, inc, Earth)
    mu = Earth.mu; J2 = Earth.J2; Re = Earth.radius;
    p = a*(1 - e^2);
    n = sqrt(mu/a^3);
    rate = 0.75*J2*n*(Re/p)^2*(5*cos(inc)^2 - 1);
end

function coe2 = coast_orbit_J2(coe1, dt, Earth)
    coe2 = coe1;
    coe2.raan = wrap2pi(coe1.raan + secular_raan_rate(coe1.a, coe1.e, coe1.i, Earth)*dt);
    coe2.argp = wrap2pi(coe1.argp + secular_argp_rate(coe1.a, coe1.e, coe1.i, Earth)*dt);

    n = sqrt(Earth.mu/coe1.a^3);
    M0 = true2mean(coe1.nu, coe1.e);
    M1 = wrap2pi(M0 + n*dt);
    coe2.nu = mean2true(M1, coe1.e);
end

function x = kep2rv_state(coe, mu)
    [r,v] = kep2rv(coe.a, coe.e, coe.i, coe.raan, coe.argp, coe.nu, mu);
    x = [r; v];
end

function mee = rv2mee(x, mu)
    coe = rv2coe(x, mu);
    mee = coe2mee(coe.a, coe.e, coe.i, coe.raan, coe.argp, coe.nu);
end

function mee = coe2mee(a, e, inc, raan, argp, nu)
    p = a*(1 - e^2);
    f = e*cos(raan + argp);
    g = e*sin(raan + argp);
    h = tan(inc/2)*cos(raan);
    k = tan(inc/2)*sin(raan);
    L = wrap2pi(raan + argp + nu);
    mee = [p; f; g; h; k; L];
end

function B = mee_B(mee, mu)
    p = mee(1); f = mee(2); g = mee(3); h = mee(4); k = mee(5); L = mee(6);

    w  = 1 + f*cos(L) + g*sin(L);
    s2 = 1 + h^2 + k^2;
    q  = sqrt(p/mu);

    B = zeros(6,3);

    B(1,:) = [0,                    2*p/w*q,                          0];
    B(2,:) = [q*sin(L),             q*((w+1)*cos(L)+f)/w,           -q*(h*sin(L)-k*cos(L))/w];
    B(3,:) = [-q*cos(L),            q*((w+1)*sin(L)+g)/w,            q*(h*cos(L)+k*sin(L))/w];
    B(4,:) = [0,                    0,                                q*s2*cos(L)/(2*w)];
    B(5,:) = [0,                    0,                                q*s2*sin(L)/(2*w)];
    B(6,:) = [0,                    0,                                q*(h*sin(L)-k*cos(L))/w];
end

function coe = rv2coe(x, mu)
    r = x(1:3);
    v = x(4:6);

    rmag = norm(r);
    vmag = norm(v);

    hvec = cross(r,v);
    hmag = norm(hvec);

    nvec = cross([0;0;1], hvec);
    nmag = norm(nvec);

    evec = ((vmag^2 - mu/rmag)*r - dot(r,v)*v)/mu;
    e = norm(evec);

    eps = vmag^2/2 - mu/rmag;
    a = -mu/(2*eps);

    i = acos(hvec(3)/hmag);

    if nmag < 1e-12
        raan = 0;
    else
        raan = atan2(nvec(2), nvec(1));
    end

    if e < 1e-12 || nmag < 1e-12
        argp = 0;
    else
        argp = atan2(dot(cross(nvec,evec), hvec)/hmag, dot(nvec,evec));
    end

    if e < 1e-12
        nu = atan2(dot(cross(nvec,r), hvec)/hmag, dot(nvec,r));
    else
        nu = atan2(dot(cross(evec,r), hvec)/hmag, dot(evec,r));
    end

    coe.a    = a;
    coe.e    = e;
    coe.i    = i;
    coe.raan = wrap2pi(raan);
    coe.argp = wrap2pi(argp);
    coe.nu   = wrap2pi(nu);
end

function [r,v] = kep2rv(a,e,inc,raan,argp,nu,mu)
    p = a*(1 - e^2);

    rpf = (p/(1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
    vpf = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];

    cO = cos(raan); sO = sin(raan);
    ci = cos(inc);  si = sin(inc);
    co = cos(argp); so = sin(argp);

    Q = [ cO*co - sO*so*ci,   -cO*so - sO*co*ci,   sO*si;
          sO*co + cO*so*ci,   -sO*so + cO*co*ci,  -cO*si;
          so*si,               co*si,              ci   ];

    r = Q*rpf;
    v = Q*vpf;
end

function M = true2mean(nu, e)
    E = 2*atan2(sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2));
    M = wrap2pi(E - e*sin(E));
end

function nu = mean2true(M, e)
    E = solve_kepler(M, e);
    nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    nu = wrap2pi(nu);
end

function E = solve_kepler(M, e)
    M = angwrap(M);

    if e < 0.8
        E = M;
    else
        E = pi;
    end

    for k = 1:50
        f  = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        dE = -f/fp;
        E  = E + dE;
        if abs(dE) < 1e-13
            break
        end
    end
end

function a = angwrap(a)
    a = mod(a + pi, 2*pi) - pi;
end

function a = wrap2pi(a)
    a = mod(a, 2*pi);
    if a < 0
        a = a + 2*pi;
    end
end