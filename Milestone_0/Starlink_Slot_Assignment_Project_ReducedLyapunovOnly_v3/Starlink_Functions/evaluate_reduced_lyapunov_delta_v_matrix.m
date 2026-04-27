function [DV, solGrid, estGrid] = evaluate_reduced_lyapunov_delta_v_matrix(cfg)
%EVALUATE_REDUCED_LYAPUNOV_DELTA_V_MATRIX Simulate every satellite-slot pair.

    nSat = numel(cfg.sats);
    nSlot = numel(cfg.slots);

    firstSol = simulate_reduced_lyapunov_slot(cfg.sats(1), cfg.slots(1), cfg);
    solGrid = repmat(firstSol, nSat, nSlot);

    firstEst = firstSol.plan;
    estGrid = repmat(firstEst, nSat, nSlot);

    total = zeros(nSat, nSlot);
    raw = zeros(nSat, nSlot);
    cleanup = zeros(nSat, nSlot);
    tf_days = zeros(nSat, nSlot);
    raanErr_deg = zeros(nSat, nSlot);
    mErr_deg = zeros(nSat, nSlot);
    success = false(nSat, nSlot);

    for i = 1:nSat
        for j = 1:nSlot
            if i == 1 && j == 1
                sol = firstSol;
            else
                fprintf('  Evaluating %s against %s with reduced Lyapunov ...\n', ...
                    cfg.sats(i).id, cfg.slots(j).id);
                sol = simulate_reduced_lyapunov_slot(cfg.sats(i), cfg.slots(j), cfg);
            end

            solGrid(i,j) = sol;
            estGrid(i,j) = sol.plan;

            total(i,j) = sol.deltaV_mps;
            raw(i,j) = sol.rawDeltaV_mps;
            cleanup(i,j) = sol.cleanupDeltaV_mps;
            tf_days(i,j) = sol.tf/86400;
            raanErr_deg(i,j) = rad2deg(sol.finalRAANError);
            mErr_deg(i,j) = rad2deg(sol.finalMError);
            success(i,j) = sol.success;
        end
    end

    DV = struct();
    DV.total_mps = total;
    DV.raw_mps = raw;
    DV.cleanup_mps = cleanup;
    DV.tf_days = tf_days;
    DV.raanErr_deg = raanErr_deg;
    DV.mErr_deg = mErr_deg;
    DV.success = success;
end
