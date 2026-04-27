function InterFamily = get_intermediary_orbit_family(Earth, Parking, Target, t_coast, N_orbits)

% get_intermediary_orbit_family
%
% Builds a family of intermediary orbits that all achieve the required RAAN
% drift over the chosen coast time.
%
% Inputs:
%   Earth   - struct with Earth.J2, Earth.radius, Earth.mu
%   Parking - parking orbit struct
%   Target  - target orbit struct
%   t_coast - coast time [s]
%
% Output:
%   InterFamily - 1xN struct array of intermediary orbit candidates

    %% User settings
    ecc_inter = 1e-3;

    % Practical SMA cutoffs
    sma_min_user = Parking.a;   % must be above parking orbit
    sma_max_user = 3.0 * Target.a;      % ceiling to avoid huge drift orbits

    %% Constants
    J2     = Earth.J2;
    rEarth = Earth.radius;
    mu     = Earth.mu;

    %% Required RAAN drift rate
    raanDot_T = SecularDrift(Earth, Target);

    deltaRAAN0 = angdiff(deg2rad(Parking.raan), deg2rad(Target.raan));

    raanDot_req = raanDot_T + deltaRAAN0 / t_coast;

    %% Feasibility bounds from |cos(i)| <= 1
    %
    % J2 RAAN drift:
    %
    %   OmegaDot = -(3/2) J2 n (Re/p)^2 cos(i)
    %
    % Solving for cos(i):
    %
    %   cos(i) = -(2/3) * OmegaDot/J2 * p^2/(n Re^2)
    %
    % With p = a(1-e^2), n = sqrt(mu/a^3),
    %
    %   cos(i) = C * a^(7/2)
    %
    % where C is constant for fixed e and required RAAN drift.

    C = -(2/3) * (raanDot_req / J2) * ...
        ((1 - ecc_inter^2)^2 / (sqrt(mu) * rEarth^2));

    if abs(C) == 0
        error('Required RAAN drift is zero. Inclination is not uniquely constrained by SMA.');
    end

    a_feas_max = (1 / abs(C))^(2/7);

    sma_lower = sma_min_user;
    sma_upper = min(a_feas_max, sma_max_user);

    if sma_lower >= sma_upper
        error(['No feasible SMA range exists with the current constraints.\n', ...
               'Try increasing sma_max_user, changing t_coast, or relaxing the parking-orbit cutoff.']);
    end

    %% Sample SMA family
    a_vec = linspace(sma_upper, sma_lower, N_orbits);

    %% Build family
    InterFamily = repmat(struct(), 1, N_orbits);

    for k = 1:N_orbits

        a = a_vec(k);
        ecc = ecc_inter;

        p = a * (1 - ecc^2);
        n = sqrt(mu / a^3);

        cosi = -(2/3) * (raanDot_req / J2) * ...
               (p^2 / (n * rEarth^2));

        % Numerical safety
        if abs(cosi) > 1
            warning('Skipping orbit %d because |cos(i)| > 1.', k);
            continue
        end

        InterFamily(k).a       = a;
        InterFamily(k).ecc     = ecc;
        InterFamily(k).inc     = acosd(cosi);
        InterFamily(k).raan    = Parking.raan;
        InterFamily(k).argp    = Parking.argp;
        InterFamily(k).M       = Parking.M;

        InterFamily(k).n       = n;
        InterFamily(k).p       = p;
        InterFamily(k).raanDot = raanDot_req;

        InterFamily(k).alt_perigee = a * (1 - ecc) - rEarth;
        InterFamily(k).alt_apogee  = a * (1 + ecc) - rEarth;

    end

    %% Print summary
    fprintf('\nIntermediary orbit family generated:\n');
    fprintf('Required RAAN drift rate: %.12e rad/s\n', raanDot_req);
    fprintf('SMA lower bound: %.3f km\n', sma_lower);
    fprintf('SMA upper bound: %.3f km\n', sma_upper);
    fprintf('\n');

    fprintf('  k        a [km]       inc [deg]     hp [km]      ha [km]\n');
    fprintf('-------------------------------------------------------------\n');

    for k = 1:N_orbits
        fprintf('%3d   %12.3f   %10.5f   %10.3f   %10.3f\n', ...
            k, ...
            InterFamily(k).a, ...
            InterFamily(k).inc, ...
            InterFamily(k).alt_perigee, ...
            InterFamily(k).alt_apogee);
    end

end