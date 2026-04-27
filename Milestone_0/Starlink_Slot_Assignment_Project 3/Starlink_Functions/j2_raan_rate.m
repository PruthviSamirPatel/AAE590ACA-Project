function raanDot = j2_raan_rate(a, inc_rad, Earth)
%J2_RAAN_RATE Secular J2 nodal regression rate for circular orbit.
%
%   raanDot = -3/2 J2 n (Re/a)^2 cos(i)
%
% Units:
%   a        [km], scalar or array
%   inc_rad  [rad]
%   raanDot  [rad/s]
%
% Vectorized: MATLAB integral() may pass a vector of semi-major axes.

    n = sqrt(Earth.mu ./ (a.^3));
    raanDot = -1.5 .* Earth.J2 .* n .* (Earth.radius ./ a).^2 .* cos(inc_rad);
end
