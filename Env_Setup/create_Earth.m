rEarth =  6378.1363; % km
mu = 398600.4415; % km3/s2

Earth.radius = rEarth;
Earth.mu = mu;
Earth.J2 = 0.001082635; 

save('Earth_params.mat', 'Earth')