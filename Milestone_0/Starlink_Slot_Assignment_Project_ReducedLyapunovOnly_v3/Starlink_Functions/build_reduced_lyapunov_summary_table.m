function T = build_reduced_lyapunov_summary_table(cfg, best, solGrid)
%BUILD_REDUCED_LYAPUNOV_SUMMARY_TABLE Summary for chosen assignment.

    n = numel(best.assignment);

    Sat = strings(n,1);
    Slot = strings(n,1);
    ArrivalEpoch = strings(n,1);
    Coast_days = zeros(n,1);
    Controlled_days = zeros(n,1);
    Tf_days = zeros(n,1);
    RawDV_mps = zeros(n,1);
    CleanupDV_mps = zeros(n,1);
    TotalDV_mps = zeros(n,1);
    TotalDV_WithMargin_mps = zeros(n,1);
    FinalAError_km = zeros(n,1);
    FinalRAANError_deg = zeros(n,1);
    FinalMResidual_deg = zeros(n,1);
    MinAltitude_km = zeros(n,1);
    Status = strings(n,1);

    for k = 1:n
        i = best.assignment(k).satIdx;
        j = best.assignment(k).slotIdx;
        sol = solGrid(i,j);

        Sat(k) = string(cfg.sats(i).id);
        Slot(k) = string(cfg.slots(j).id);
        ArrivalEpoch(k) = string(char(sol.arrivalEpoch));
        Coast_days(k) = sol.tCoast/86400;
        Controlled_days(k) = sol.tControlled/86400;
        Tf_days(k) = sol.tf/86400;
        RawDV_mps(k) = sol.rawDeltaV_mps;
        CleanupDV_mps(k) = sol.cleanupDeltaV_mps;
        TotalDV_mps(k) = sol.deltaV_mps;
        TotalDV_WithMargin_mps(k) = best.assignment(k).deltaV_withMargin_mps;
        FinalAError_km(k) = sol.finalAError_km;
        FinalRAANError_deg(k) = rad2deg(sol.finalRAANError);
        FinalMResidual_deg(k) = rad2deg(sol.finalMError);
        MinAltitude_km(k) = sol.minAltitude_km;
        Status(k) = string(char(sol.status));
    end

    T = table(Sat, Slot, ArrivalEpoch, Coast_days, Controlled_days, Tf_days, ...
        RawDV_mps, CleanupDV_mps, TotalDV_mps, TotalDV_WithMargin_mps, ...
        FinalAError_km, FinalRAANError_deg, FinalMResidual_deg, MinAltitude_km, Status);
end
