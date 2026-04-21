function Xdot = EoM_minTime_barrier(~, X, B, mu, umax, rSafe, rOn, rhoPenalty, epsBar)

    % =========================
    % Unpack State
    % =========================
    r = X(1:3);        % position 
    v = X(4:6);        % velocity 
    lamb = X(7:12);    % costate

    % =========================
    % Optimal Control
    % =========================
    p = -B' * lamb;                % switching function

    if norm(p) < 1e-10
        u = zeros(3,1);            % avoid numerical issues
    else
        u = (p / norm(p)) * umax;
    end

    % =========================
    % State Dynamics
    % =========================
    r_norm = norm(r);

    a_grav = -mu * r / r_norm^3;

    f0 = [v;
          a_grav];

    xDot = f0 + B*u;

    % =========================
    % Costate Dynamics
    % =========================
    G = 3*mu*(r*r')/r_norm^5 - mu*eye(3)/r_norm^3;

    A = [zeros(3), eye(3);
         G,        zeros(3)];

    % running cost penalty gradient
    [~, gradPhi] = radiusPenaltyCost(r, rSafe, rOn, epsBar);

    lambDot = -A' * lamb - [rhoPenalty * gradPhi;
                            zeros(3,1)];

    % =========================
    % Output
    % =========================
    Xdot = [xDot; lambDot];

end