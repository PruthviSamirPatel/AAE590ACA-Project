clc; clear; close all

%% ===================== USER SETTINGS =====================
days_coast = 10;
t_coast = days_coast*24*60*60; % s

N_intermediary = 50;

tFinal_leg1_days = 10;
tFinal_leg2_days = 5;

rTol_km   = 1.0;      % rendezvous position tolerance
vTol_km_s = 1e-3;     % rendezvous velocity tolerance

saveOnlyConverged = true;

%% ===================== LOAD ENVIRONMENT =====================
load("Earth_params.mat")
rEarth = Earth.radius;
mu     = Earth.mu;

load("Orbital_Slot.mat")
load("Parking_Orbit.mat")

%% ===================== NONDIMENSIONALIZATION =====================
lStar = rEarth;
tStar = sqrt(rEarth^3/mu);
vStar = lStar/tStar;
aStar = lStar/tStar^2;

mu_nd = 1;
t_coast_nd = t_coast/tStar;

rTol_nd = rTol_km/lStar;
vTol_nd = vTol_km_s/vStar;

%% ===================== CONTROL SETTINGS =====================
uMax_km_s2 = 1e-3;
uMax = uMax_km_s2/aStar;

%% ===================== EARTH AVOIDANCE SETTINGS =====================
hSafe_km = 150;
hOn_km   = 300;

rSafe = (rEarth + hSafe_km)/lStar;
rOn   = (rEarth + hOn_km)/lStar;

kBarrier = 5e-5;
epsBar   = 1e-6;

%% ===================== GAINS =====================
Kmee = eye(5);
Pmee = eye(5);

Kcart = eye(3);
Pcart = eye(3);

%% ===================== CONSTANTS =====================
Re_nd = Earth.radius/lStar;
J2 = Earth.J2;

%% ===================== BUILD INTERMEDIARY FAMILY =====================
InterFamily = get_intermediary_orbit_family(Earth, Parking, Target, ...
                                            t_coast, N_intermediary);

N = length(InterFamily);

ConvergedTransfers = struct([]);
AllResults = {};

fprintf('\n============================================================\n');
fprintf('Beginning intermediary-orbit batch solve with %d candidates\n', N);
fprintf('============================================================\n');

iGood = 0;

%% ===================== LOOP OVER INTERMEDIARY ORBITS =====================
for kk = 1:N

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Solving candidate %d of %d\n', kk, N);
    fprintf('a = %.3f km, i = %.6f deg\n', InterFamily(kk).a, InterFamily(kk).inc);
    fprintf('------------------------------------------------------------\n');

    Result = struct( ...
    'caseID', kk, ...
    'InterDesired', InterFamily(kk), ...
    'success', false, ...
    'leg1Success', false, ...
    'leg2Success', false, ...
    'failureReason', "", ...
    'DV1_m_s', NaN, ...
    'DV2_m_s', NaN, ...
    'DV_total_m_s', NaN, ...
    'final_dr_km', NaN, ...
    'final_dv_km_s', NaN, ...
    'minAlt_leg1_km', NaN, ...
    'minAlt_leg2_km', NaN ...
    );

    try
        %% ================================================================
        %% LEG 1: PARKING TO INTERMEDIARY ORBIT USING MEE LYAPUNOV CONTROL
        %% ================================================================

        Inter = InterFamily(kk);

        % Initial orbit
        a0    = Parking.a;
        e0    = Parking.ecc;
        inc0  = deg2rad(Parking.inc);
        raan0 = deg2rad(Parking.raan);
        argp0 = deg2rad(Parking.argp);
        M0    = deg2rad(Parking.M);

        % Intermediary orbit
        aI    = Inter.a;
        eI    = Inter.ecc;
        incI  = deg2rad(Inter.inc);
        raanI = deg2rad(Inter.raan);
        argpI = deg2rad(Inter.argp);

        % Initial true anomaly
        E0 = kepler(rad2deg(M0), e0);
        ta0 = 2*atan2(sqrt(1+e0)*sind(E0/2), ...
                      sqrt(1-e0)*cosd(E0/2));

        % Target true anomaly for MEE slow-element targeting
        taI = 0;

        % Convert to MEE
        mee0_dim = kep2mee(a0, e0, inc0, raan0, argp0, ta0);
        meeT_dim = kep2mee(aI, eI, incI, raanI, argpI, taI);

        X0_leg1 = [mee0_dim(1)/lStar; mee0_dim(2:6)];
        xslowT_nd = [meeT_dim(1)/lStar; meeT_dim(2:5)];

        tSpan_leg1 = [0, tFinal_leg1_days*24*3600/tStar];

        opts1 = odeset('RelTol',1e-12, 'AbsTol',1e-12, ...
            'Events', @(t,X) earth_keepout_event_mee(t, X, rSafe));

        [t1, X1, te1, Xe1, ie1] = ode45(@(t,X) mee_lyap_dyn(t, X, xslowT_nd, ...
            Kmee, Pmee, mu_nd, uMax, rSafe, rOn, kBarrier, epsBar, J2, Re_nd), ...
            tSpan_leg1, X0_leg1, opts1);

        % Check if Leg 1 hit Earth
        if ~isempty(ie1)
            Result.failureReason = "Leg 1 crossed Earth keepout.";
            Result.leg1EventID = ie1(end);
            % AllResults = [AllResults; Result];
            AllResults{end+1,1} = Result;
            fprintf('FAILED: %s\n', Result.failureReason);
            continue
        end

        % Recover achieved intermediary orbit
        mee_dim_end = [X1(end,1)*lStar; X1(end,2:6)'];
        [a_end, e_end, inc_end, raan_end, argp_end, ta_end] = mee2kep(mee_dim_end);

        InterAchieved.a    = a_end;
        InterAchieved.ecc  = e_end;
        InterAchieved.inc  = rad2deg(inc_end);
        InterAchieved.raan = rad2deg(raan_end);
        InterAchieved.argp = rad2deg(argp_end);

        E_end = 2*atan2(sqrt(1-e_end)*sin(ta_end/2), ...
                        sqrt(1+e_end)*cos(ta_end/2));
        E_end = mod(E_end, 2*pi);

        M_end = E_end - e_end*sin(E_end);
        M_end = mod(M_end, 2*pi);

        InterAchieved.M = rad2deg(M_end);
        InterAchieved.n = sqrt(Earth.mu/InterAchieved.a^3);

        % Leg 1 control and Delta-V
        num1 = length(t1);
        u1_hist = zeros(num1,3);
        u1_norm = zeros(num1,1);
        alt1_hist = zeros(num1,1);

        for ii = 1:num1
            [~, u_out] = mee_lyap_dyn_plot_helper(t1(ii), X1(ii,:)', ...
                xslowT_nd, Kmee, Pmee, mu_nd, uMax, rSafe, rOn, ...
                kBarrier, epsBar, J2, Re_nd);

            u1_hist(ii,:) = u_out'*aStar;
            u1_norm(ii) = norm(u1_hist(ii,:));

            p_i = X1(ii,1);
            f_i = X1(ii,2);
            g_i = X1(ii,3);
            L_i = X1(ii,6);
            w_i = 1 + f_i*cos(L_i) + g_i*sin(L_i);
            r_i = p_i/w_i*lStar;
            alt1_hist(ii) = r_i - rEarth;
        end

        DV1_km_s = trapz(t1*tStar, u1_norm);
        DV1_m_s = 1000*DV1_km_s;

        % if min(alt1_hist) < hSafe_km
        %     Result.failureReason = "Leg 1 altitude went below safe altitude.";
        %     AllResults = [AllResults; Result];
        %     fprintf('FAILED: %s\n', Result.failureReason);
        %     continue
        % end

        Result.leg1Success = true;

        %% ================================================================
        %% COAST TO RAAN MATCH, THEN FIND CLOSEST PHASING POINT
        %% ================================================================
        
        [tMatch, raanI_match, raanT_match] = find_raan_match_time(Earth, ...
            InterAchieved, Target);
        
        fprintf('RAAN match coast = %.6f days\n', tMatch/(24*3600));
        
        % Search extra coast after RAAN match
        N_extra_revs = 5;       % coast this many target revolutions after RAAN match
        N_phase_grid = 500;     % grid resolution before fminbnd refinement
        
        TargetAtMatch = propagate_orbit_J2(Earth, Target, tMatch);
        Ttarget = 2*pi/sqrt(Earth.mu/TargetAtMatch.a^3);
        
        dtExtraMax = N_extra_revs*Ttarget;
        
        [dtExtraClosest, minSep_km, minRelSpeed_km_s] = find_closest_extra_coast( ...
            Earth, InterAchieved, Target, tMatch, dtExtraMax, N_phase_grid);
        
        tCoastTotal = tMatch + dtExtraClosest;
        
        InterCoasted  = propagate_orbit_J2(Earth, InterAchieved, tCoastTotal);
        TargetCoasted = propagate_orbit_J2(Earth, Target,        tCoastTotal);
        
        fprintf('Extra coast to closest approach = %.6f days\n', dtExtraClosest/(24*3600));
        fprintf('Total pre-Leg-2 coast           = %.6f days\n', tCoastTotal/(24*3600));
        fprintf('Initial Leg 2 separation        = %.6f km\n', minSep_km);
        fprintf('Initial Leg 2 relative speed    = %.9f km/s\n', minRelSpeed_km_s);
        fprintf('RAAN error at Leg 2 start       = %.9f deg\n', ...
            wrapTo180Local(InterCoasted.raan - TargetCoasted.raan));

        %% ================================================================
        %% LEG 2: INTERMEDIARY TO TARGET SLOT USING CARTESIAN LYAPUNOV
        %% ================================================================

        % Chaser initial orbit after coast
        a0    = InterCoasted.a;
        e0    = InterCoasted.ecc;
        inc0  = InterCoasted.inc;
        raan0 = InterCoasted.raan;
        argp0 = InterCoasted.argp;
        M0    = InterCoasted.M;

        % Target after same coast
        aT    = TargetCoasted.a;
        eT    = TargetCoasted.ecc;
        incT  = TargetCoasted.inc;
        raanT = TargetCoasted.raan;
        argpT = TargetCoasted.argp;
        MT    = TargetCoasted.M;

        % Convert mean anomaly to true anomaly
        E0_deg = kepler(M0, e0);
        ta0_deg = 2*atan2d(sqrt(1+e0)*sind(E0_deg/2), ...
                           sqrt(1-e0)*cosd(E0_deg/2));
        ta0_deg = mod(ta0_deg, 360);

        ET_deg = kepler(MT, eT);
        taT_deg = 2*atan2d(sqrt(1+eT)*sind(ET_deg/2), ...
                           sqrt(1-eT)*cosd(ET_deg/2));
        taT_deg = mod(taT_deg, 360);

        % Cartesian states
        [x,y,z,vx,vy,vz] = kep2cart(a0, e0, inc0, argp0, raan0, ta0_deg, mu);
        rC0_nd = [x;y;z]/lStar;
        vC0_nd = [vx;vy;vz]/vStar;

        [x,y,z,vx,vy,vz] = kep2cart(aT, eT, incT, argpT, raanT, taT_deg, mu);
        rT0_nd = [x;y;z]/lStar;
        vT0_nd = [vx;vy;vz]/vStar;

        X0_leg2 = [rC0_nd; vC0_nd; rT0_nd; vT0_nd];

        tSpan_leg2 = [0, tFinal_leg2_days*24*3600/tStar];

        opts2 = odeset('RelTol',1e-12, 'AbsTol',1e-12, ...
            'Events', @(t,X) combined_events_cart(t, X, rSafe, rTol_nd, vTol_nd));

        [t2, X2, te2, Xe2, ie2] = ode45(@(t,X) cart_lyap_dyn(t, X, ...
            Kcart, Pcart, mu_nd, uMax, rSafe, rOn, kBarrier, epsBar, J2, Re_nd), ...
            tSpan_leg2, X0_leg2, opts2);

        rC_nd = X2(:,1:3);
        vC_nd = X2(:,4:6);
        rT_nd = X2(:,7:9);
        vT_nd = X2(:,10:12);

        dr_km   = (rC_nd - rT_nd)*lStar;
        dv_km_s = (vC_nd - vT_nd)*vStar;

        final_dr_km = norm(dr_km(end,:));
        final_dv_km_s = norm(dv_km_s(end,:));

        num2 = length(t2);
        u2_hist = zeros(num2,3);
        u2_norm = zeros(num2,1);
        alt2_hist = zeros(num2,1);

        for ii = 1:num2
            r_i = X2(ii,1:3)';
            alt2_hist(ii) = norm(r_i)*lStar - rEarth;

            [~, u_out] = cart_lyap_dyn_plot_helper(t2(ii), X2(ii,:)', ...
                Kcart, Pcart, mu_nd, uMax, rSafe, rOn, kBarrier, epsBar, J2, Re_nd);

            u2_hist(ii,:) = u_out'*aStar;
            u2_norm(ii) = norm(u2_hist(ii,:));
        end

        DV2_km_s = trapz(t2*tStar, u2_norm);
        DV2_m_s = 1000*DV2_km_s;

        % if min(alt2_hist) < hSafe_km
        %     Result.failureReason = "Leg 2 altitude went below safe altitude.";
        %     AllResults = [AllResults; Result];
        %     fprintf('FAILED: %s\n', Result.failureReason);
        %     continue
        % end

        leg2HitEarth = false;
        leg2Rendezvous = false;

        if ~isempty(ie2)
            if ie2(end) == 1
                leg2HitEarth = true;
            elseif ie2(end) == 2
                leg2Rendezvous = true;
            end
        end

        if leg2HitEarth
            Result.failureReason = "Leg 2 crossed Earth keepout.";
            AllResults = [AllResults; Result];
            fprintf('FAILED: %s\n', Result.failureReason);
            continue
        end

        if final_dr_km <= rTol_km && final_dv_km_s <= vTol_km_s
            leg2Rendezvous = true;
        end

        if ~leg2Rendezvous
            Result.failureReason = "Leg 2 did not converge to rendezvous tolerance.";
            Result.final_dr_km = final_dr_km;
            Result.final_dv_km_s = final_dv_km_s;
            AllResults = [AllResults; Result];

            fprintf('FAILED: %s\n', Result.failureReason);
            fprintf('Final ||dr|| = %.6f km\n', final_dr_km);
            fprintf('Final ||dv|| = %.9f km/s\n', final_dv_km_s);
            continue
        end

        Result.leg2Success = true;
        Result.success = true;
        Result.failureReason = "None";
        Result.InterAchieved = InterAchieved;
        Result.InterCoasted = InterCoasted;
        Result.TargetCoasted = TargetCoasted;
        
        Result.tMatch_s = tMatch;
        Result.tMatch_days = tMatch/(24*3600);
        Result.raanI_match_deg = raanI_match;
        Result.raanT_match_deg = raanT_match;
        %% ================================================================
        %% STORE SUCCESSFUL TRANSFER
        %% ================================================================

        Result.InterAchieved = InterAchieved;
        Result.InterCoasted = InterCoasted;
        Result.TargetCoasted = TargetCoasted;

        Result.tMatch_s = tMatch;
        Result.tMatch_days = tMatch/(24*3600);
        Result.raanI_match_deg = raanI_match;
        Result.raanT_match_deg = raanT_match;

        Result.DV1_m_s = DV1_m_s;
        Result.DV2_m_s = DV2_m_s;
        Result.DV_total_m_s = DV1_m_s + DV2_m_s;

        Result.final_dr_km = final_dr_km;
        Result.final_dv_km_s = final_dv_km_s;

        Result.minAlt_leg1_km = min(alt1_hist);
        Result.minAlt_leg2_km = min(alt2_hist);

        Result.t1_days = t1*tStar/(24*3600);
        Result.X1 = X1;
        Result.u1_hist_km_s2 = u1_hist;
        Result.alt1_hist_km = alt1_hist;

        Result.t2_days = t2*tStar/(24*3600);
        Result.X2 = X2;
        Result.u2_hist_km_s2 = u2_hist;
        Result.alt2_hist_km = alt2_hist;

        Result.rC_km = X2(:,1:3)*lStar;
        Result.vC_km_s = X2(:,4:6)*vStar;
        Result.rT_km = X2(:,7:9)*lStar;
        Result.vT_km_s = X2(:,10:12)*vStar;
        Result.dr_km = dr_km;
        Result.dv_km_s = dv_km_s;

        AllResults = [AllResults; Result];

        iGood = iGood + 1;
        ConvergedTransfers(iGood).caseID = kk;
        ConvergedTransfers(iGood).InterDesired = Result.InterDesired;
        ConvergedTransfers(iGood).InterAchieved = Result.InterAchieved;
        ConvergedTransfers(iGood).InterCoasted = Result.InterCoasted;
        ConvergedTransfers(iGood).TargetCoasted = Result.TargetCoasted;

        ConvergedTransfers(iGood).tMatch_days = Result.tMatch_days;
        ConvergedTransfers(iGood).DV1_m_s = Result.DV1_m_s;
        ConvergedTransfers(iGood).DV2_m_s = Result.DV2_m_s;
        ConvergedTransfers(iGood).DV_total_m_s = Result.DV_total_m_s;

        ConvergedTransfers(iGood).final_dr_km = Result.final_dr_km;
        ConvergedTransfers(iGood).final_dv_km_s = Result.final_dv_km_s;
        ConvergedTransfers(iGood).minAlt_leg1_km = Result.minAlt_leg1_km;
        ConvergedTransfers(iGood).minAlt_leg2_km = Result.minAlt_leg2_km;

        ConvergedTransfers(iGood).t1_days = Result.t1_days;
        ConvergedTransfers(iGood).X1 = Result.X1;
        ConvergedTransfers(iGood).t2_days = Result.t2_days;
        ConvergedTransfers(iGood).X2 = Result.X2;
        ConvergedTransfers(iGood).rC_km = Result.rC_km;
        ConvergedTransfers(iGood).rT_km = Result.rT_km;
        ConvergedTransfers(iGood).dr_km = Result.dr_km;
        ConvergedTransfers(iGood).dv_km_s = Result.dv_km_s;

        fprintf('SUCCESS.\n');
        fprintf('Leg 1 DV = %.6f m/s\n', DV1_m_s);
        fprintf('Leg 2 DV = %.6f m/s\n', DV2_m_s);
        fprintf('Total DV = %.6f m/s\n', DV1_m_s + DV2_m_s);
        fprintf('Final ||dr|| = %.6f km\n', final_dr_km);
        fprintf('Final ||dv|| = %.9f km/s\n', final_dv_km_s);

    catch ME
        Result.success = false;
        Result.failureReason = string(ME.message);
        AllResults = [AllResults; Result];

        fprintf('FAILED WITH ERROR:\n');
        fprintf('%s\n', ME.message);
        continue
    end

end

%% ===================== POSTPROCESS SUCCESSFUL CASES =====================
fprintf('\n============================================================\n');
fprintf('Batch solve complete.\n');
fprintf('Number converged: %d out of %d\n', length(ConvergedTransfers), N);
fprintf('============================================================\n');

if ~isempty(ConvergedTransfers)

    DV_all = [ConvergedTransfers.DV_total_m_s];
    [DV_best, idxBest] = min(DV_all);

    BestTransfer = ConvergedTransfers(idxBest);

    fprintf('\nBest converged transfer:\n');
    fprintf('Case ID = %d\n', BestTransfer.caseID);
    fprintf('Total DV = %.6f m/s\n', DV_best);
    fprintf('Leg 1 DV = %.6f m/s\n', BestTransfer.DV1_m_s);
    fprintf('Leg 2 DV = %.6f m/s\n', BestTransfer.DV2_m_s);
    fprintf('Intermediary SMA = %.6f km\n', BestTransfer.InterDesired.a);
    fprintf('Intermediary inc = %.6f deg\n', BestTransfer.InterDesired.inc);
    fprintf('Final ||dr|| = %.6f km\n', BestTransfer.final_dr_km);
    fprintf('Final ||dv|| = %.9f km/s\n', BestTransfer.final_dv_km_s);

else
    BestTransfer = [];
    fprintf('\nNo intermediary orbit candidates converged.\n');
end

save('Batch_Intermediary_Transfer_Results.mat', ...
     'InterFamily', 'AllResults', 'ConvergedTransfers', 'BestTransfer');

%% ===================== SUMMARY PLOT =====================
if ~isempty(ConvergedTransfers)

    caseIDs = [ConvergedTransfers.caseID];
    sma_all = arrayfun(@(s) s.InterDesired.a, ConvergedTransfers);
    inc_all = arrayfun(@(s) s.InterDesired.inc, ConvergedTransfers);
    DV_all = [ConvergedTransfers.DV_total_m_s];

    figure('Color','w');
    plot(sma_all, DV_all, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor','k');
    grid on;
    xlabel('Intermediary SMA [km]');
    ylabel('Total \DeltaV [m/s]');
    title('Converged Transfers: Total \DeltaV vs Intermediary SMA');

    figure('Color','w');
    plot(inc_all, DV_all, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor','k');
    grid on;
    xlabel('Intermediary Inclination [deg]');
    ylabel('Total \DeltaV [m/s]');
    title('Converged Transfers: Total \DeltaV vs Intermediary Inclination');

end

function [dtBest, minSep_km, minRelSpeed_km_s] = find_closest_extra_coast( ...
    Earth, InterAchieved, Target, tMatch, dtExtraMax, Ngrid)

    dtGrid = linspace(0, dtExtraMax, Ngrid);
    sepGrid = zeros(size(dtGrid));

    for i = 1:length(dtGrid)
        [sepGrid(i), ~] = separation_after_extra_coast( ...
            dtGrid(i), tMatch, Earth, InterAchieved, Target);
    end

    [~, idxBest] = min(sepGrid);

    % Build local refinement interval around best grid point
    idxA = max(idxBest - 1, 1);
    idxB = min(idxBest + 1, length(dtGrid));

    dtA = dtGrid(idxA);
    dtB = dtGrid(idxB);

    if dtA == dtB
        dtBest = dtGrid(idxBest);
    else
        costFun = @(dt) separation_after_extra_coast( ...
            dt, tMatch, Earth, InterAchieved, Target);

        dtBest = fminbnd(costFun, dtA, dtB);
    end

    [minSep_km, minRelSpeed_km_s] = separation_after_extra_coast( ...
        dtBest, tMatch, Earth, InterAchieved, Target);
end

function [sep_km, relSpeed_km_s] = separation_after_extra_coast( ...
    dtExtra, tMatch, Earth, InterAchieved, Target)

    tNow = tMatch + dtExtra;

    InterNow  = propagate_orbit_J2(Earth, InterAchieved, tNow);
    TargetNow = propagate_orbit_J2(Earth, Target,        tNow);

    [rI, vI] = orbit_struct_to_cart(InterNow, Earth.mu);
    [rT, vT] = orbit_struct_to_cart(TargetNow, Earth.mu);

    sep_km = norm(rI - rT);
    relSpeed_km_s = norm(vI - vT);
end

function [r, v] = orbit_struct_to_cart(Orb, mu)

    a    = Orb.a;
    e    = Orb.ecc;
    inc  = Orb.inc;
    raan = Orb.raan;
    argp = Orb.argp;
    M    = Orb.M;

    E_deg = kepler(M, e);

    ta_deg = 2*atan2d(sqrt(1+e)*sind(E_deg/2), ...
                      sqrt(1-e)*cosd(E_deg/2));

    ta_deg = mod(ta_deg, 360);

    [x, y, z, vx, vy, vz] = kep2cart(a, e, inc, argp, raan, ta_deg, mu);

    r = [x; y; z];
    v = [vx; vy; vz];
end

function ang = wrapTo180Local(ang)
    ang = mod(ang + 180, 360) - 180;
end