function [x, y, z, vx, vy, vz] = kep2cart(a, ecc, inc, argp, raan, ta, mu)
% This program converts from keplerian orbital elements to cartesian elements
% 
% INPUTS:
% a: semi-major axis [km]
% ecc: eccentricity
% inc: inclination [deg]
% argp: argument of periapsis [deg]
% raan: right ascension of ascending node [deg]
% ta: true anomaly [deg]
% mu: gravitational parameter [km3/s2]
%
% OUTPUT: position (x,y,z) and velocity (vx,vy,vz) in ECI frame [km] or
% [km/s]
% 
% Pruthvi Samir Patel

p = a*(1-ecc^2); % semilatus rectum
% Rotation matrix:
R = [[cosd(raan)*cosd(argp) - sind(raan)*sind(argp)*cosd(inc), ...
-cosd(raan)*sind(argp) - sind(raan)*cosd(argp)*cosd(inc), sind(raan)*sind(inc)]; ...
[sind(raan)*cosd(argp) + cosd(raan)*sind(argp)*cosd(inc), ...
-sind(raan)*sind(argp) + cosd(raan)*cosd(argp)*cosd(inc), -cosd(raan)*sind(inc)]; ...
[sind(argp)*sind(inc), cosd(argp)*sind(inc), cosd(inc)]];
% Position vector in perifocal frame
r = p/(1+ecc*cosd(ta)); 
r_peri = r*[cosd(ta); sind(ta); 0]; % [e^, theta^, h^]
% Velocity vector in perifocal frame
v = sqrt(mu/(a*(1-ecc^2)));
v_peri = v*[-sind(ta); ecc+cosd(ta); 0];
% Convert frames
r_cart = R*r_peri;
v_cart = R*v_peri;
% Assign to return variables:
x = r_cart(1);
y = r_cart(2);
z = r_cart(3);
vx = v_cart(1);
vy = v_cart(2);
vz = v_cart(3);
end
