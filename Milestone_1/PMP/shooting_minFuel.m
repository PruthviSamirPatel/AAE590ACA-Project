%% FCTN
function Psi = shooting_minFuel(lambda0, X0, t0, tf, umax, mu_nd, B, xTargetf, rho, opts)
    % setup initial state
    X_init = [X0;lambda0];
    tSpan = [t0, tf];

    % propagate
    traj = ode45(@(t,X) EoM_minFuel(t, X, B, mu_nd, umax, rho), tSpan, X_init, opts);

    % get final state:
    X_his = deval(traj, tf);
    xf = X_his(1:6);
    Psi = xf - xTargetf;
end