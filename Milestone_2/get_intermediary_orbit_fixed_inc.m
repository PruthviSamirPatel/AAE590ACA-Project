function Inter = get_intermediary_orbit_fixed_inc(Earth, Parking, Target, t_coast)
% operates in dimensional form
% Chooses an intermediary orbit with fixed inclination and small e
% so that RAAN matches the target after t_coast.
%
% Solves for SMA using J2 secular RAAN drift:
%
% OmegaDot = -(3/2)*J2*n*(Re/p)^2*cos(i)
%
% where p = a*(1 - e^2), n = sqrt(mu/a^3)

    J2     = Earth.J2;
    rEarth = Earth.radius;
    mu     = Earth.mu;

    %% Choose intermediary orbit properties
    Inter.ecc  = 1e-3;
    Inter.inc  = Target.inc;      % choose this fixed inclination [deg]
    Inter.raan = Parking.raan;
    Inter.argp = Parking.argp;
    Inter.M    = Parking.M;

    %% Required RAAN drift rate
    raanDot_T = SecularDrift(Earth, Target);

    deltaRAAN0 = angdiff(deg2rad(Parking.raan), deg2rad(Target.raan));

    raanDot_req = raanDot_T + deltaRAAN0/t_coast;

    %% Solve for SMA
    inc_rad = deg2rad(Inter.inc);
    e = Inter.ecc;

    % From:
    % raanDot_req = -(3/2)*J2*sqrt(mu/a^3)*(rEarth/(a*(1-e^2)))^2*cos(i)
    %
    % Therefore:
    % raanDot_req = C * a^(-7/2)
    %
    % where:
    % C = -(3/2)*J2*sqrt(mu)*rEarth^2*cos(i)/(1-e^2)^2

    C = -(3/2)*J2*sqrt(mu)*rEarth^2*cos(inc_rad)/(1 - e^2)^2;

    if sign(C) ~= sign(raanDot_req)
        error('No feasible SMA exists: chosen inclination gives RAAN drift with wrong sign.');
    end

    a_sol = (C/raanDot_req)^(2/7);

    if ~isreal(a_sol) || a_sol <= rEarth
        error('Solved SMA is nonphysical. Check chosen inclination or coast time.');
    end

    %% Store result
    Inter.a = a_sol;
    Inter.n = sqrt(mu/Inter.a^3);
    Inter.raanDot = raanDot_req;

    fprintf('\nIntermediary orbit solved with fixed inclination:\n');
    fprintf('  a        = %.6f km\n', Inter.a);
    fprintf('  altitude = %.6f km\n', Inter.a - rEarth);
    fprintf('  ecc      = %.6f\n', Inter.ecc);
    fprintf('  inc      = %.6f deg\n', Inter.inc);
    fprintf('  RAANdot  = %.12e rad/s\n', Inter.raanDot);
end