function r = circular_elements_to_rv(a, inc, raan, M)
%CIRCULAR_ELEMENTS_TO_RV Position vector for circular orbit.
%
% Inputs:
%   a    [km]
%   inc  [rad]
%   raan [rad]
%   M    [rad], equal to true anomaly for circular orbit.
%
% Argument of perigee is arbitrary for e=0 and set to zero.

    u = M;

    r_pf = [a*cos(u); a*sin(u); 0];

    cO = cos(raan); sO = sin(raan);
    ci = cos(inc);  si = sin(inc);

    R3O = [ cO, -sO, 0;
            sO,  cO, 0;
             0,   0, 1];

    R1i = [1,  0,   0;
           0, ci, -si;
           0, si,  ci];

    r = R3O * R1i * r_pf;
end
