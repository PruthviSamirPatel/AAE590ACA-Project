function [a, ecc, inc, argp, raan, ta] = cart2kep(x, y, z, vx, vy, vz, mu)
% This program converts from cartesian elements (ECI)
% to keplerian or orbital elements
% 
% INPUTS: 
% position (x,y,z) in ECI frame [km]
% velocity (vx,vy,vz) in ECI frame [km/s]
% mu: gravitational parameter [km3/s2]
%
% OUTPUTS:
% a: semi-major axis [km]
% ecc: eccentricity
% inc: inclination [deg]
% argp: argument of periapsis [deg]
% raan: right ascension of ascending node [deg]
% ta: true anomaly [deg]
%
% Pruthvi Samir Patel

% Position and velocity
r_cart = [x; y; z]; % position in cartesian coordinates
v_cart = [vx; vy; vz]; % velocity in cartesian coordinates
r = norm(r_cart);
% Specific angular momentum
h_cart = cross(r_cart, v_cart); 
h = norm(h_cart);
hHat = h_cart/h; % unit vector
% Eccentricity vector
ecc_cart = cross(v_cart,h_cart)/mu - r_cart/r;
ecc = norm(ecc_cart);
eccHat = ecc_cart/ecc;
% Semimajor axis
a = h^2/(mu*(1-ecc^2)); 
% Inclination
inc = acosd(hHat(3)); 
% Line of nodes
nodes = [-hHat(2); hHat(1); 0];
nodesHat = nodes/norm(nodes);
% Right ascension of ascending node
raan = atan2d(hHat(1),-hHat(2));
% Argument of periapsis
cosw = dot(nodesHat,eccHat);
sinw = dot(cross(nodesHat,eccHat), hHat); 
argp = atan2d(sinw,cosw); 
% True anomaly
cosf = dot(eccHat, r_cart/r);
sinf = dot(cross(eccHat,r_cart/r), hHat); 
ta = atan2d(sinf,cosf);
end
