%% EoM FUNCTION
function Xdot = EoM_minFuel(~, X, B, mu, umax, rho)
    % breakup state
    r = X(1:3); v = X(4:6);
    lamb = X(7:12);

    % get optimal control
    p = -B'*lamb;
    p_norm = norm(p);
    uHat = p / p_norm;
    
    % Switching Function: S > 0 when ||p|| > 1
    S = p_norm - 1; 
    
    Gamma = umax/2 * (1 + tanh(S/rho));
    u = uHat * Gamma;

    % State dynamics
    r_norm = norm(r);
    a_grav = -mu * r / r_norm^3;
    f0 = [v;
          a_grav];
    xDot = f0 + B*u;

    % Costate Dynamics
    G = 3*mu*(r*r')/r_norm^5 - mu*eye(3)/r_norm^3;
    A = [zeros(3), eye(3);
         G,        zeros(3)];
    lambDot = -A' * lamb;

    Xdot = [xDot; lambDot];
end
