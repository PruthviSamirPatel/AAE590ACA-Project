function [value, isterminal, direction] = event_target_or_earth(t, X, P)

r = X(1:3);
v = X(4:6);

[rT, vT, ~, ~] = targetStateAtTime(t, P.Target, P.mu_nd);

e_r = norm(r - rT);
e_v = norm(v - vT);

% Event 1: close enough to target trajectory/state
tol_r = 5e-3;
tol_v = 5e-3;
value1 = max(e_r - tol_r, e_v - tol_v);

% Event 2: Earth collision
value2 = norm(r) - P.r_min;

value = [value1; value2];
isterminal = [1; 1];
direction = [-1; -1];
end