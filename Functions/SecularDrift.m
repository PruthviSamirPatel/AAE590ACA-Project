%% Secular J2 RAAN drift
function raanDot = SecularDrift(Earth, Orbit)
    J2 = Earth.J2;
    rEarth = Earth.radius;
    n = Orbit.n;
    p = Orbit.a * (1 - Orbit.ecc^2);
    inc = Orbit.inc;

    raanDot = -3/2 * J2 * (rEarth^2 / p^2) * n * cosd(inc);
end