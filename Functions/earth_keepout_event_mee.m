%% Event function
function [value, isterminal, direction] = earth_keepout_event_mee(~, X, rSafe)
    p = X(1);
    f = X(2);
    g = X(3);
    L = X(6);

    w = 1 + f*cos(L) + g*sin(L);

    if p <= 0 || w <= 0
        value = -1;
        isterminal = 1;
        direction = -1;
        return
    end

    r = p / w;

    value = r - rSafe;
    isterminal = 1;
    direction = -1;
end