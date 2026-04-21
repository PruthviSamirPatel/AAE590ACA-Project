function Psi = shooting_minTime_barrier(Z, x0, t0, umax, mu_nd, B, Target, ...
                                        opts, rSafe, rOn, rhoPenalty, epsBar)

lambda0 = Z(1:6);
tf = Z(7);

if tf <= t0
    Psi = 1e6 * ones(7,1);
    return
end

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
xDotTargetf = [vxTf; vyTf; vzTf; -mu_nd*[xTf; yTf; zTf]/norm([xTf; yTf; zTf])^3];

%% Propagate Spacecraft
X0 = [x0; lambda0];

try
    sol = ode45(@(t,X) EoM_minTime_barrier(t, X, B, mu_nd, umax, ...
        rSafe, rOn, rhoPenalty, epsBar), [t0 tf], X0, opts);

    % reject if integration terminated early
    if sol.x(end) < tf - 1e-10
        Psi = 1e6 * ones(7,1);
        return
    end

    Xf = deval(sol, tf);

catch
    Psi = 1e6 * ones(7,1);
    return
end

xf = Xf(1:6);
lf = Xf(7:12);

%% Hamiltonian at tf
rf = xf(1:3); 
vf = xf(4:6);

p_f = -B'*lf; 
if norm(p_f) < 1e-10
    uf = zeros(3,1);
else
    uf = (p_f / norm(p_f)) * umax;
end

f0_f = [vf; -mu_nd*rf/norm(rf)^3];
[phi_f, ~] = radiusPenaltyCost(rf, rSafe, rOn, epsBar);

Hf = 1 + rhoPenalty*phi_f + lf' * (f0_f + B*uf);

%% Transversality Conditions
Psi = [xf - xTargetf;
       Hf - lf' * xDotTargetf];
end