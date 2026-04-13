function Psi = shooting_minTime(Z, x0, t0, umax, mu_nd, B, Target, opts)
lambda0 = Z(1:6);
tf = Z(7);

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
sol = ode45(@(t,X) EoM_minTime(t, X, B, mu_nd, umax), [t0 tf], X0, opts);
Xf = deval(sol, tf);

xf = Xf(1:6);
lf = Xf(7:12);

%% Hamiltonian at tf
rf = xf(1:3); vf = xf(4:6);
p_f = -B'*lf; 
if norm(p_f) < 1e-10
    uf = zeros(3,1);
else
    uf = (p_f / norm(p_f)) * umax;
end
f0_f = [vf; -mu_nd*rf/norm(rf)^3];
Hf = 1 + lf' * (f0_f + B*uf);

%% Transversality Conditions
Psi = [xf - xTargetf;                % Boundary condition: hit target orbit state
       Hf - lf' * xDotTargetf];      % Transversality: H(tf) - lambda'*xDot_target = 0
end