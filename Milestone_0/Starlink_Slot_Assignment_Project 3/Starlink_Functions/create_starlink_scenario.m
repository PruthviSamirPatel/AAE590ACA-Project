function cfg = create_starlink_scenario()
%CREATE_STARLINK_SCENARIO  Starlink-like 5-satellite / 5-slot scenario.
%
% This scenario is intentionally small and deterministic.
%
% Starlink-like assumptions:
%   - low Earth orbit shell,
%   - fixed inclination,
%   - low-thrust tangential orbit raising,
%   - J2 differential nodal precession for RAAN-plane deployment.
%
% Only SMA and RAAN are changed by control/planning. Mean anomaly is used
% for slot assignment and residual reporting.

    cfg = struct();

    %% Earth
    cfg.Earth.mu     = 398600.4418;        % km^3/s^2
    cfg.Earth.radius = 6378.1363;          % km
    cfg.Earth.J2     = 1.08262668e-3;      % -

    %% Epochs
    cfg.launchEpoch = datetime(2026,5,1,0,0,0,'TimeZone','UTC');

    %% Starlink-like shell geometry
    cfg.hParking_km = 350.0;               % staging / drift shell altitude [km]
    cfg.hFinal_km   = 550.0;               % operational shell altitude [km]
    cfg.inc_deg     = 53.05;               % representative Starlink shell inclination [deg]
    cfg.inc_rad     = deg2rad(cfg.inc_deg);

    cfg.aParking = cfg.Earth.radius + cfg.hParking_km;
    cfg.aFinal   = cfg.Earth.radius + cfg.hFinal_km;

    %% Control
    % 2e-7 km/s^2 = 2e-4 m/s^2 = 0.2 mm/s^2.
    % This produces a multi-day low-thrust raise, similar in structure to
    % electric orbit raising. Change this value for your vehicle model.
    cfg.uMax_km_s2 = 2.0e-7;

    %% Constraints
    cfg.aMin = cfg.aParking;               % do not dip below drift/parking SMA
    cfg.aMax = cfg.aFinal + 80.0;          % optional guard for Lyapunov overshoot
    cfg.aTol_km = 0.5;
    cfg.raanTol_deg = 0.02;

    %% Assignment weights
    cfg.assignment.timeWeight = 1.0;        % days
    cfg.assignment.meanAnomalyWeight_daysPerDeg = 0.015;
    cfg.assignment.raanErrorWeight_daysPerDeg = 100.0;

    %% Lyapunov settings
    cfg.lyap.kA = 1.0e-5;                  % feedback gain in km/s^2 per km before saturation
    cfg.lyap.extraTimeFraction = 0.20;     % max simulation time beyond PMP tf
    cfg.lyap.maxExtraDays = 5.0;
    cfg.lyap.dtPlot = 600.0;               % seconds
    cfg.lyap.tolA_km = 0.5;
    cfg.lyap.tolRAAN_deg = 0.05;

    %% Initial five satellites
    launchRAAN_deg = 40.0;
    satMeanAnom_deg = [0, 72, 144, 216, 288];

    for i = 1:5
        cfg.sats(i).id = sprintf('SAT%02d', i);
        cfg.sats(i).a0 = cfg.aParking;
        cfg.sats(i).raan0_deg = launchRAAN_deg;
        cfg.sats(i).M0_deg = satMeanAnom_deg(i);
    end

    %% Design five open slots.
    % Two slots share the same operational RAAN plane.
    % Three are lower in RAAN, reached by waiting longer at the lower SMA.
    %
    % The "same RAAN" plane is defined as the plane reached by immediate
    % max-thrust raising from parking to final shell. This is physically
    % consistent because even an immediate low-thrust raise accumulates
    % some differential J2 precession relative to a satellite already at
    % the final shell.

    [~, dOmRaise, dMRaise] = compute_raise_integrals(cfg.aParking, cfg.aFinal, cfg.uMax_km_s2, cfg);
    baseOffset_deg = rad2deg(dOmRaise);    % immediate-raise operational RAAN offset

    lowerStep_deg = 1.0;
    slotOffsets_deg = [baseOffset_deg, ...
                       baseOffset_deg, ...
                       baseOffset_deg - lowerStep_deg, ...
                       baseOffset_deg - 2*lowerStep_deg, ...
                       baseOffset_deg - 3*lowerStep_deg];

    % Choose open slot mean anomalies in a realistic map. These are
    % generated from the nominal assigned satellites plus small offsets so
    % the assignment problem is not trivial but remains well-conditioned.
    nominalSatForSlot = [1, 2, 3, 4, 5];
    meanPerturb_deg = [8, -10, 6, -7, 12];

    raanDotP = j2_raan_rate(cfg.aParking, cfg.inc_rad, cfg.Earth);
    raanDotF = j2_raan_rate(cfg.aFinal,   cfg.inc_rad, cfg.Earth);
    nP = circ_mean_rate(cfg.aParking, cfg.Earth);
    nF = circ_mean_rate(cfg.aFinal,   cfg.Earth);

    for j = 1:5
        slotRAAN0_deg = launchRAAN_deg + slotOffsets_deg(j);
        D0 = wrapToPiLocal(deg2rad(launchRAAN_deg - slotRAAN0_deg));

        tCoast = -(D0 + dOmRaise)/(raanDotP - raanDotF);
        if tCoast < 0
            tCoast = 0;
        end

        satIdx = nominalSatForSlot(j);
        Mslot0_deg = cfg.sats(satIdx).M0_deg + rad2deg((nP-nF)*tCoast + dMRaise) + meanPerturb_deg(j);

        cfg.slots(j).id = sprintf('SLOT%02d', j);
        cfg.slots(j).aF = cfg.aFinal;
        cfg.slots(j).raan0_deg = mod(slotRAAN0_deg, 360);
        cfg.slots(j).M0_deg = mod(Mslot0_deg, 360);
        cfg.slots(j).nominalCoast_s = tCoast;
        cfg.slots(j).nominalEpoch = cfg.launchEpoch + seconds(tCoast + compute_raise_time(cfg.aParking, cfg.aFinal, cfg.uMax_km_s2, cfg));
    end
end
