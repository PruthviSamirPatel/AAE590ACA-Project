%% load environment
load('Earth_params.mat');
rEarth = Earth.radius;
mu = Earth.mu;

%% Make target orbit at initial epoch
alt = 630; 

Target.a = rEarth + alt; % semi-major axis
Target.ecc = 0.0005; % eccentricity
Target.inc = 51.9; % inclination [deg]
Target.raan = 50; % right ascension of ascending node [deg]
Target.argp = 0; % argument of periapsis [deg]
Target.M = 145; % mean anomaly [deg]

save('Orbital_Slot.mat', 'Target')