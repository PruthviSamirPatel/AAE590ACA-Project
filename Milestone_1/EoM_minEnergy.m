function Xdot = EoM_minEnergy(~, X, B, mu)
    % breakup state
    x = X(1:4);
    r = x(1:2); v = x(3:4);
    lamb = X(5:end);

    % get optimal control
    u = -1/2 * B' * lamb;

    % get derivative of state
    f0 = [v; -mu*r/norm(r)^3];
    xDot = f0 + B*u;

    % get derivative of costate
    G = 3*mu*r*r'/norm(r)^5 - mu*eye(2)/norm(r)^3; % gravity gradient wrt pos
    A = [zeros(2), eye(2);
         G, zeros(2)]; % gradient of f0 wrt full state x
    lambDot = -A' * lamb;

    % return state derivative
    Xdot = [xDot; lambDot];
end