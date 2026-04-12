rEarth =  6378.1363; % km
mu = 398600.4415; % km3/s2

Earth.radius = rEarth;
Earth.mu = mu;

save('Earth_params.mat', 'Earth')