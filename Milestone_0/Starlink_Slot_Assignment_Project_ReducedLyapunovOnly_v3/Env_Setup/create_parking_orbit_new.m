%% load environment
load('Earth_params.mat');
rEarth = Earth.radius;
mu = Earth.mu;

%% Make parking orbit at initial epoch
alt = 610; 

Parking.a = rEarth + alt; % semi-major axis
Parking.ecc = 0.001; % eccentricity
Parking.inc = 51.9; % inclination [deg]
Parking.raan = 60; % right ascension of ascending node [deg]
Parking.argp = 0; % argument of periapsis [deg]
Parking.M = 10; % mean anomaly [deg]

save('Parking_Orbit_new.mat', 'Parking')