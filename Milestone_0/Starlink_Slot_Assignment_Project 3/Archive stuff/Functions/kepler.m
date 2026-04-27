%% kepler
function [E] = kepler(M,ecc)
% this function estimates the Eccentric Anomaly [deg] value given Mean Anomaly
% [deg] through an iterative process following the Netwon-Rhapsod
% method 
% made by Pruthvi Patel
tol = 1e-12;
% to start, I will let E=M
M = deg2rad(M);
Enew = M;
Eold = 0;
while (abs(Enew-Eold)>tol)
    Eold = Enew;
    Enew = Eold - (Eold - ecc*sin(Eold)-M)/(1-ecc*cos(Eold));
end
E = rad2deg(Enew);
% f = 2*atand(sqrt((1+ecc)/(1-ecc))*tand(E/2));
end