function [value, isterminal, direction] = earth_keepout_event_OC(~, X, rSafe)

r = X(1:3);

value = norm(r) - rSafe;
isterminal = 1;
direction = -1;

end