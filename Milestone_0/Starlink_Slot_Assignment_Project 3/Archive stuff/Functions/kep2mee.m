%% Keplerian -> MEE
function mee = kep2mee(a,e,i,raan,argp,ta)
    p = a*(1 - e^2);

    f = e*cos(raan + argp);
    g = e*sin(raan + argp);

    h = tan(i/2)*cos(raan);
    k = tan(i/2)*sin(raan);

    L = raan + argp + ta;
    L = mod(L, 2*pi);

    mee = [p; f; g; h; k; L];
end