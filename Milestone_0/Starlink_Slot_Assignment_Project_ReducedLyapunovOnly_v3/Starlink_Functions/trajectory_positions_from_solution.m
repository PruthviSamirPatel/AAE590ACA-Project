function R = trajectory_positions_from_solution(sol, sat, slot, cfg)
%TRAJECTORY_POSITIONS_FROM_SOLUTION Convert reduced element history to ECI positions.

    n = numel(sol.t);
    R = zeros(3,n);

    nF = circ_mean_rate(slot.aF, cfg.Earth);
    OmDotF = j2_raan_rate(slot.aF, cfg.inc_rad, cfg.Earth);

    for k = 1:n
        tk = sol.t(k);
        slotRAAN = deg2rad(slot.raan0_deg) + OmDotF*tk;
        slotM    = deg2rad(slot.M0_deg)    + nF*tk;

        satRAAN = slotRAAN + sol.relRAAN(k);
        satM    = slotM    + sol.relM(k);

        R(:,k) = circular_elements_to_rv(sol.a(k), cfg.inc_rad, satRAAN, satM);
    end
end
