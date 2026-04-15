function Xdot = EoM_LyapunovPenalty(t, X, P)

r = X(1:3);
v = X(4:6);

[u, ~] = LyapunovControl(t, X, P);

r_norm = norm(r);
a_grav = -P.mu_nd * r / r_norm^3;

Xdot = [v;
        a_grav + u];
end