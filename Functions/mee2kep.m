%% MEE -> Keplerian
function [a,e,i,raan,argp,ta] = mee2kep(mee)
    p = mee(1);
    f = mee(2);
    g = mee(3);
    h = mee(4);
    k = mee(5);
    L = mee(6);

    e = sqrt(f^2 + g^2);
    a = p / (1 - e^2);

    raan = atan2(k,h);
    i = 2*atan(sqrt(h^2 + k^2));

    lonper = atan2(g,f);    % Omega + omega
    argp = lonper - raan;
    ta = L - lonper;

    raan = mod(raan, 2*pi);
    argp = mod(argp, 2*pi);
    ta   = mod(ta,   2*pi);
end