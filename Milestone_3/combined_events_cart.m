function [value, isterminal, direction] = combined_events_cart(~, X, rSafe, rTol, vTol)
    % Event 1: Earth keep-out
    rC = X(1:3);
    rmag = norm(rC);

    value1 = rmag - rSafe;

    % Event 2: Rendezvous achieved
    dr = X(1:3) - X(7:9);
    dv = X(4:6) - X(10:12);

    value2 = max(norm(dr) - rTol, norm(dv) - vTol);

    value = [value1; value2];
    isterminal = [1; 1];
    direction  = [-1; -1];
end