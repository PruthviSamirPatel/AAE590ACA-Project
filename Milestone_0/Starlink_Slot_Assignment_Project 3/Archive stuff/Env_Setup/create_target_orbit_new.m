%% load environment
load('Earth_params.mat');
rEarth = Earth.radius;
mu = Earth.mu;

%% Make target orbit at initial epoch
alt = 630; 

a = rEarth + alt; % semi-major axis
ecc = 0.0005; % eccentricity
inc = 51.9; % inclination [deg]
raan = 40; % right ascension of ascending node [deg]
argp = 0; % argument of periapsis [deg]
M = 30; % mean anomaly [deg]

n = sqrt(mu/a^3); % mean motion

Target.a = a;
Target.ecc = ecc;
Target.inc = inc;
Target.raan = raan;
Target.argp = argp;
Target.M = M; 
Target.n = n; 

%% Nondimensionalize:
% Get nondim quantities
lStar = rEarth;
tStar = sqrt(rEarth^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;
mu_nd = 1;

% add nonDim:
Target.a_nd = a/lStar;
Target.n_nd = n*tStar;

save('Orbital_Slot_new.mat', 'Target')