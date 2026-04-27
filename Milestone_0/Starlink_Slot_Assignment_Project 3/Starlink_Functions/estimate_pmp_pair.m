function est = estimate_pmp_pair(sat, slot, cfg)
%ESTIMATE_PMP_PAIR Fast reduced PMP estimate for one satellite-slot pair.
%
% Structure:
%   - Coast at aParking to collect differential J2 RAAN drift.
%   - Raise from aParking to aFinal with maximum tangential acceleration.
%
% This is the reduced minimum-time path-constrained PMP structure when:
%   1) a >= aParking,
%   2) lower RAAN is desired,
%   3) RAAN drift is fastest at the lower permitted SMA.

    [tRaise, dOmRaise, dMRaise] = compute_raise_integrals(sat.a0, slot.aF, cfg.uMax_km_s2, cfg);

    raanDot0 = j2_raan_rate(sat.a0, cfg.inc_rad, cfg.Earth);
    raanDotF = j2_raan_rate(slot.aF, cfg.inc_rad, cfg.Earth);

    n0 = circ_mean_rate(sat.a0, cfg.Earth);
    nF = circ_mean_rate(slot.aF, cfg.Earth);

    D0 = wrapToPiLocal(deg2rad(sat.raan0_deg - slot.raan0_deg));
    M0rel = wrapToPiLocal(deg2rad(sat.M0_deg - slot.M0_deg));

    denom = raanDot0 - raanDotF;

    if abs(denom) < 1e-14
        tCoast = 0;
        status = "No differential RAAN drift";
    else
        tCoast = -(D0 + dOmRaise)/denom;
        status = "OK";
        if tCoast < 0
            % Negative coast time means the slot is not lower enough for
            % this drift-only strategy. The fastest admissible plan is
            % immediate raise, with nonzero RAAN residual.
            tCoast = 0;
            status = "Immediate raise; RAAN residual remains";
        end
    end

    tf = tCoast + tRaise;

    finalRAANError = wrapToPiLocal(D0 + (raanDot0 - raanDotF)*tCoast + dOmRaise);
    finalMError    = wrapToPiLocal(M0rel + (n0 - nF)*tCoast + dMRaise);

    est.tCoast = tCoast;
    est.tRaise = tRaise;
    est.tf = tf;
    est.finalRAANError = finalRAANError;
    est.finalMError = finalMError;
    est.deltaV_mps = 1000 * cfg.uMax_km_s2 * tRaise;
    est.status = status;
end
