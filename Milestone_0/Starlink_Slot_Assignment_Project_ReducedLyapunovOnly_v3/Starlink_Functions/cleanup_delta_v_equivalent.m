function cleanupDV = cleanup_delta_v_equivalent(raanErr, mErr, cfg)
%CLEANUP_DELTA_V_EQUIVALENT Engineering equivalent DV for residual slot error.
%
% This is a reporting/assignment metric. It is not an impulsive maneuver model.

    cleanupDV = cfg.cleanup.raan_mpsPerDeg * abs(rad2deg(wrapToPiLocal(raanErr))) + ...
                cfg.cleanup.meanAnomaly_mpsPerDeg * abs(rad2deg(wrapToPiLocal(mErr)));
end
