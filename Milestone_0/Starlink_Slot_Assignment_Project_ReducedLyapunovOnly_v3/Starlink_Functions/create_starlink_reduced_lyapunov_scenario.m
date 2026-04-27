function cfg = create_starlink_reduced_lyapunov_scenario()
%CREATE_STARLINK_REDUCED_LYAPUNOV_SCENARIO Deterministic 5-sat/5-slot case.
%
% Reduced Lyapunov-only version. There are no PMP settings.
%
% IMPORTANT:
%   uMax_km_s2 is an acceleration, not a velocity.  The default is the
%   Project-3 low-thrust value 0.2 mm/s^2.  A value of 1 m/s^2 would make
%   the 200 km altitude raise take only about 112 seconds, which appears
%   vertical on multi-day plots.

    cfg = struct();

    %% Earth
    cfg.Earth.mu     = 398600.4418;        % km^3/s^2
    cfg.Earth.radius = 6378.1363;          % km
    cfg.Earth.J2     = 1.08262668e-3;      % -

    %% Epoch
    cfg.launchEpoch = datetime(2026,5,1,0,0,0,'TimeZone','UTC');

    %% Shell geometry
    cfg.hParking_km = 350.0;
    cfg.hFinal_km   = 550.0;
    cfg.inc_deg     = 53.05;
    cfg.inc_rad     = deg2rad(cfg.inc_deg);

    cfg.aParking = cfg.Earth.radius + cfg.hParking_km;
    cfg.aFinal   = cfg.Earth.radius + cfg.hFinal_km;

    %% Control: Project-3 low-thrust default
    cfg.uMax_km_s2 = 2.0e-7;               % 0.0002 m/s^2 = 0.2 mm/s^2

    %% Bounds / tolerances
    cfg.aMin = cfg.aParking;
    cfg.aMax = cfg.aFinal + 25.0;
    cfg.aTol_km = 0.5;
    cfg.raanTol_deg = 0.05;

    %% Lyapunov settings
    cfg.lyap.kA = 1.0e-5;                  % km/s^2 per km before saturation
    cfg.lyap.tolA_km = 0.5;
    cfg.lyap.extraRaiseTimeFraction = 0.30;
    cfg.lyap.maxExtraDays = 3.0;
    cfg.lyap.dtPlot = 600.0;               % output grid seconds

    %% Cleanup equivalent Delta-V model for residual slot errors
    cfg.cleanup.meanAnomaly_mpsPerDeg = 0.05; % 8 deg -> 0.4 m/s
    cfg.cleanup.raan_mpsPerDeg = 3.0;

    %% Assignment
    cfg.assignment.dvMarginFraction = 0.20;

    %% Satellites
    launchRAAN_deg = 40.0;
    satMeanAnom_deg = [0, 72, 144, 216, 288];

    for i = 1:5
        cfg.sats(i).id = sprintf('SAT%02d', i);
        cfg.sats(i).a0 = cfg.aParking;
        cfg.sats(i).raan0_deg = launchRAAN_deg;
        cfg.sats(i).M0_deg = satMeanAnom_deg(i);
    end

    %% Slots: two in same RAAN, three lower RAAN.
    % The first RAAN plane corresponds to immediate low-thrust raise from
    % parking to final shell.  Three additional slots are 1, 2, and 3 deg
    % lower in RAAN and are reached by natural J2 drift during the coast.
    [~, dOmRaise, dMRaise] = compute_raise_integrals(cfg.aParking, cfg.aFinal, cfg.uMax_km_s2, cfg);

    baseOffset_deg = rad2deg(dOmRaise);
    lowerStep_deg = 1.0;
    slotOffsets_deg = [baseOffset_deg, ...
                       baseOffset_deg, ...
                       baseOffset_deg - lowerStep_deg, ...
                       baseOffset_deg - 2*lowerStep_deg, ...
                       baseOffset_deg - 3*lowerStep_deg];

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
        Mslot0_deg = cfg.sats(satIdx).M0_deg + ...
            rad2deg((nP-nF)*tCoast + dMRaise) + meanPerturb_deg(j);

        cfg.slots(j).id = sprintf('SLOT%02d', j);
        cfg.slots(j).aF = cfg.aFinal;
        cfg.slots(j).raan0_deg = mod(slotRAAN0_deg, 360);
        cfg.slots(j).M0_deg = mod(Mslot0_deg, 360);
        cfg.slots(j).nominalCoast_s = tCoast;
        cfg.slots(j).nominalEpoch = cfg.launchEpoch + ...
            seconds(tCoast + compute_raise_time(cfg.aParking, cfg.aFinal, cfg.uMax_km_s2, cfg));
    end
end
