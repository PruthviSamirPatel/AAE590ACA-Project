function est = estimate_reduced_lyapunov_pair(sat, slot, cfg)
%ESTIMATE_REDUCED_LYAPUNOV_PAIR Plan coast time for RAAN matching.
%
% This is not PMP.  It computes a coast time such that
%
%   DeltaOmega0 + (OmegaDotP - OmegaDotF)tCoast + DeltaOmegaRaise ~= 0.
%
% The actual closed-loop propagation is performed by
% simulate_reduced_lyapunov_slot.

    D0 = wrapToPiLocal(deg2rad(sat.raan0_deg - slot.raan0_deg));

    [tRaise, dOmRaise, dMRaise] = compute_raise_integrals(sat.a0, slot.aF, cfg.uMax_km_s2, cfg);

    rateP = j2_raan_rate(sat.a0, cfg.inc_rad, cfg.Earth);
    rateF = j2_raan_rate(slot.aF, cfg.inc_rad, cfg.Earth);
    rateDiff = rateP - rateF;

    if abs(rateDiff) < 1e-14
        tCoast = 0;
        status = "No differential RAAN drift";
    else
        tCoast = -(D0 + dOmRaise)/rateDiff;
        status = "OK";
        if tCoast < 0
            tCoast = 0;
            status = "Immediate raise; RAAN residual remains";
        end
    end

    nP = circ_mean_rate(sat.a0, cfg.Earth);
    nF = circ_mean_rate(slot.aF, cfg.Earth);

    M0 = wrapToPiLocal(deg2rad(sat.M0_deg - slot.M0_deg));
    OmFinal = wrapToPiLocal(D0 + rateDiff*tCoast + dOmRaise);
    MFinal = wrapToPiLocal(M0 + (nP-nF)*tCoast + dMRaise);

    rawDV = 1000*abs(cfg.uMax_km_s2)*tRaise;
    cleanupDV = cleanup_delta_v_equivalent(OmFinal, MFinal, cfg);

    est = struct();
    est.satId = sat.id;
    est.slotId = slot.id;
    est.tCoast = tCoast;
    est.tRaise = tRaise;
    est.tf = tCoast + tRaise;
    est.rawDeltaV_mps = rawDV;
    est.cleanupDeltaV_mps = cleanupDV;
    est.deltaV_mps = rawDV + cleanupDV;
    est.finalRAANError = OmFinal;
    est.finalMError = MFinal;
    est.status = status;
end
