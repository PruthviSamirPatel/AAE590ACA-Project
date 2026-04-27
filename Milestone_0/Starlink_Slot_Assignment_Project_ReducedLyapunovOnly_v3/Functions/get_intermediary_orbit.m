function Inter = get_intermediary_orbit(Earth, Parking, Target, t_coast)
% operates in dimensional form
% Chooses an intermediary orbit with fixed a and small e, then solves for i
% so that RAAN matches the target after t_coast.

%% unpack
J2 = Earth.J2;
rEarth = Earth.radius;
mu = Earth.mu;

%% chosen intermediary shape
Inter.ecc = 1e-2;                 % nearly circular
Inter.a   = 1.1 * Target.a;        % user-chosen trial SMA
Inter.raan = Parking.raan;         % starts with parking RAAN
Inter.argp = Parking.argp;

%% intermediary quantities
p = Inter.a * (1 - Inter.ecc^2);
n = sqrt(mu / Inter.a^3);

%% target secular drift
raanDot_T = SecularDrift(Earth, Target);   % rad/s

%% required intermediary drift
deltaRAAN0 = angdiff(deg2rad(Parking.raan), deg2rad(Target.raan)); % rad
raanDot_req = raanDot_T + deltaRAAN0 / t_coast;                    % rad/s

%% solve for inclination
cosi = -(2/3) * (raanDot_req / J2) * (p^2 / (n * rEarth^2));

if abs(cosi) > 1
    error('No feasible inclination exists for this chosen SMA and coast time.');
end

Inter.inc = acosd(cosi);
Inter.n = n;
Inter.raanDot = raanDot_req;
end


function raanDot = SecularDrift(Earth, Orbit)
% operates in dimensional form

J2 = Earth.J2;
rEarth = Earth.radius;
n = Orbit.n;
p = Orbit.a * (1 - Orbit.ecc^2);
inc = Orbit.inc;

raanDot = -3/2 * J2 * (rEarth^2 / p^2) * n * cosd(inc);   % rad/s
end