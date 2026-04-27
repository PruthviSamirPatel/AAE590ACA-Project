function [tRaise, dOmegaRel, dMRel] = compute_raise_integrals(a0, aF, uTangential, cfg)
%COMPUTE_RAISE_INTEGRALS Differential RAAN and M accumulated during raise.
%
% The returned dOmegaRel and dMRel are relative to a reference slot already
% moving on the final shell:
%
%   dOmegaRel = integral_0^tRaise [OmegaDot(a(t)) - OmegaDot(aF)] dt
%   dMRel     = integral_0^tRaise [n(a(t))        - n(aF)]        dt
%
% This is what matters for slot matching.

    u = abs(uTangential);
    if u <= 0
        error('compute_raise_integrals: uTangential must be positive.');
    end

    tRaise = compute_raise_time(a0, aF, u, cfg);

    aLo = min(a0, aF);
    aHi = max(a0, aF);

    nF = circ_mean_rate(aF, cfg.Earth);
    OmF = j2_raan_rate(aF, cfg.inc_rad, cfg.Earth);

    % da/dt = 2*u/n(a), so dt/da = n(a)/(2*u)
    integrandOm = @(a) (j2_raan_rate(a, cfg.inc_rad, cfg.Earth) - OmF) .* ...
                       (circ_mean_rate(a, cfg.Earth) ./ (2*u));

    integrandM  = @(a) (circ_mean_rate(a, cfg.Earth) - nF) .* ...
                       (circ_mean_rate(a, cfg.Earth) ./ (2*u));

    dOmegaRel = integral(integrandOm, aLo, aHi, 'RelTol',1e-11, 'AbsTol',1e-13);
    dMRel     = integral(integrandM,  aLo, aHi, 'RelTol',1e-11, 'AbsTol',1e-13);

    if aF < a0
        dOmegaRel = -dOmegaRel;
        dMRel = -dMRel;
    end
end
