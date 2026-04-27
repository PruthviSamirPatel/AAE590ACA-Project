function M_deg = true2mean(ta_deg, ecc)

% Guard against invalid or non-elliptic eccentricity
if ~isreal(ecc) || ~isfinite(ecc)
    error('true2mean: invalid eccentricity (non-real or non-finite).');
end

if ecc < 0
    error('true2mean: eccentricity is negative.');
end

if ecc >= 1
    error('true2mean: eccentricity >= 1, orbit is not elliptical.');
end

f = deg2rad(ta_deg);

E = 2 * atan2( sqrt(max(0,1 - ecc)) * sin(f/2), ...
               sqrt(1 + ecc) * cos(f/2) );

M = E - ecc * sin(E);
M_deg = mod(rad2deg(M), 360);

end