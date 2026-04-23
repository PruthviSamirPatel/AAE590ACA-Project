function Inter = get_intermediary_orbit(Earth, Parking, Target, t_coast)
% operates in dimensional form
% Chooses an intermediary orbit with fixed a and small e
% so that RAAN matches the target after t_coast.

    J2 = Earth.J2;
    rEarth = Earth.radius;
    mu = Earth.mu;

    Inter.ecc = 1e-3;
    Inter.a   = 1.01 * Target.a;
    Inter.raan = Parking.raan;
    Inter.argp = Parking.argp;

    p = Inter.a * (1 - Inter.ecc^2);
    n = sqrt(mu / Inter.a^3);

    raanDot_T = SecularDrift(Earth, Target);

    deltaRAAN0 = angdiff(deg2rad(Parking.raan), deg2rad(Target.raan));
    raanDot_req = raanDot_T + deltaRAAN0 / t_coast;

    cosi = -(2/3) * (raanDot_req / J2) * (p^2 / (n * rEarth^2));

    if abs(cosi) > 1
        error('No feasible inclination exists for this chosen SMA and coast time.');
    end

    Inter.inc = acosd(cosi);
    Inter.n = n;
    Inter.raanDot = raanDot_req;
end