function [tMatch, raanIntFinal, raanTarFinal] = find_raan_match_time(Earth, Inter, Target)
% Solve for coast time such that propagated intermediary RAAN matches
% propagated target RAAN under secular J2 effects.
%
% Returns:
%   tMatch [s]

    % Drift rates
    raanDot_I = SecularDrift(Earth, Inter);   % rad/s
    raanDot_T = SecularDrift(Earth, Target);  % rad/s

    % Initial wrapped difference: Target - Inter
    dRAAN0 = angdiff(deg2rad(Inter.raan), deg2rad(Target.raan)); % rad

    % Need: Omega_I + raanDot_I*t = Omega_T + raanDot_T*t
    % => (raanDot_I - raanDot_T)*t = dRAAN0
    denom = raanDot_I - raanDot_T;

    if abs(denom) < 1e-14
        error('Intermediary and target RAAN drift rates are nearly equal; no unique match time.');
    end

    tMatch = dRAAN0 / denom;

    if tMatch < 0
        warning('Computed RAAN match time is negative. The first future match may require adding 2*pi / |denom|.');
        tMatch = tMatch + 2*pi/abs(denom);
    end

    raanIntFinal = mod(Inter.raan  + rad2deg(raanDot_I*tMatch), 360);
    raanTarFinal = mod(Target.raan + rad2deg(raanDot_T*tMatch), 360);
end