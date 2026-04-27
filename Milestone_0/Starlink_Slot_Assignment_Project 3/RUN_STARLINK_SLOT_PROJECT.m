clc; clear; close all;

%% ============================================================
% STARLINK-LIKE SLOT ASSIGNMENT PROJECT
%
% Five satellites are launched into the same low parking/drift shell.
% Five open slots exist in a RAAN / mean-anomaly map.
%
% Assumption requested by user:
%   Only SMA and RAAN are controlled/changed.
%
% Model:
%   x = [a; Omega; M]
%   fixed eccentricity, fixed inclination, fixed argument of perigee.
%
% Methods:
%   1) Slot assignment by brute-force 5-by-5 assignment.
%   2) Reduced path-constrained PMP:
%         coast at low SMA for RAAN precession,
%         then max tangential thrust to final SMA.
%   3) Lyapunov SMA controller following the same RAAN-drift logic.
%% ============================================================

thisFile = mfilename('fullpath');
projectRoot = fileparts(thisFile);
addpath(genpath(projectRoot));

cfg = create_starlink_scenario();

fprintf('\n============================================================\n');
fprintf(' STARLINK-LIKE SLOT ASSIGNMENT + CONTROL PROJECT\n');
fprintf('============================================================\n');
fprintf('Launch epoch                  : %s\n', char(cfg.launchEpoch));
fprintf('Parking altitude              : %.3f km\n', cfg.hParking_km);
fprintf('Final shell altitude          : %.3f km\n', cfg.hFinal_km);
fprintf('Inclination                   : %.4f deg\n', cfg.inc_deg);
fprintf('Tangential acceleration limit : %.3e km/s^2 = %.4f mm/s^2\n', ...
    cfg.uMax_km_s2, 1e6*cfg.uMax_km_s2);
fprintf('Number of satellites          : %d\n', numel(cfg.sats));
fprintf('Number of slots               : %d\n', numel(cfg.slots));

rateP = j2_raan_rate(cfg.aParking, cfg.inc_rad, cfg.Earth);
rateF = j2_raan_rate(cfg.aFinal,   cfg.inc_rad, cfg.Earth);
fprintf('\nJ2 RAAN drift at parking      : %.6f deg/day\n', rad2deg(rateP)*86400);
fprintf('J2 RAAN drift at final shell  : %.6f deg/day\n', rad2deg(rateF)*86400);
fprintf('Differential RAAN drift       : %.6f deg/day\n', rad2deg(rateP-rateF)*86400);

[Tr, dOmR, dMR] = compute_raise_integrals(cfg.aParking, cfg.aFinal, cfg.uMax_km_s2, cfg);
fprintf('\nNominal max-thrust raise time : %.6f days\n', Tr/86400);
fprintf('RAAN gained during raise relative to final shell : %.6f deg\n', rad2deg(dOmR));
fprintf('Mean anomaly gained during raise relative to final shell : %.6f deg\n', mod(rad2deg(dMR),360));

%% Print initial maps
fprintf('\n================ INITIAL SATELLITES ================\n');
for i = 1:numel(cfg.sats)
    s = cfg.sats(i);
    fprintf('%s: a0 = %.3f km, RAAN0 = %.6f deg, M0 = %.6f deg\n', ...
        s.id, s.a0, s.raan0_deg, s.M0_deg);
end

fprintf('\n================ OPEN SLOTS ================\n');
for j = 1:numel(cfg.slots)
    q = cfg.slots(j);
    fprintf('%s: aF = %.3f km, RAAN0 = %.6f deg, M0 = %.6f deg, nominal epoch = %s\n', ...
        q.id, q.aF, q.raan0_deg, q.M0_deg, char(q.nominalEpoch));
end

%% Slot assignment
[assignment, C, pairEst] = assign_slots_bruteforce(cfg.sats, cfg.slots, cfg);

fprintf('\n================ ASSIGNMENT COST MATRIX ================\n');
fprintf('Rows = satellites, columns = slots. Cost units are equivalent days.\n');
fprintf('%12s', '');
for j = 1:numel(cfg.slots)
    fprintf('%14s', cfg.slots(j).id);
end
fprintf('\n');
for i = 1:numel(cfg.sats)
    fprintf('%12s', cfg.sats(i).id);
    for j = 1:numel(cfg.slots)
        fprintf('%14.4f', C(i,j));
    end
    fprintf('\n');
end

fprintf('\n================ SELECTED SLOT ASSIGNMENT ================\n');
fprintf('%8s %8s %12s %12s %12s %12s %12s\n', ...
    'Sat','Slot','tf[d]','coast[d]','raise[d]','RAANerr[deg]','Merr[deg]');
for k = 1:numel(assignment)
    i = assignment(k).satIdx;
    j = assignment(k).slotIdx;
    e = pairEst(i,j);
    fprintf('%8s %8s %12.4f %12.4f %12.4f %12.6f %12.6f\n', ...
        cfg.sats(i).id, cfg.slots(j).id, e.tf/86400, e.tCoast/86400, ...
        e.tRaise/86400, rad2deg(e.finalRAANError), rad2deg(e.finalMError));
end

%% Run PMP and Lyapunov for assigned pairs
pmpResults = repmat(struct(), numel(assignment), 1);
lyapResults = repmat(struct(), numel(assignment), 1);

fprintf('\n================ RUNNING CONTROLLERS ================\n');

for k = 1:numel(assignment)
    i = assignment(k).satIdx;
    j = assignment(k).slotIdx;

    sat  = cfg.sats(i);
    slot = cfg.slots(j);

    fprintf('\n--- Pair %d: %s -> %s ---\n', k, sat.id, slot.id);

    pmp = solve_pmp_reduced_slot(sat, slot, cfg);
    pmpResults(k).satId = sat.id;
    pmpResults(k).slotId = slot.id;
    pmpResults(k).solution = pmp;

    fprintf('PMP reduced solution:\n');
    fprintf('  coast time       = %.6f days\n', pmp.tCoast/86400);
    fprintf('  raise time       = %.6f days\n', pmp.tRaise/86400);
    fprintf('  final time       = %.6f days\n', pmp.tf/86400);
    fprintf('  arrival epoch    = %s\n', char(pmp.arrivalEpoch));
    fprintf('  delta-V          = %.6f m/s\n', pmp.deltaV_mps);
    fprintf('  final a error    = %.6e km\n', pmp.finalAError_km);
    fprintf('  final RAAN error = %.6e deg\n', rad2deg(pmp.finalRAANError));
    fprintf('  final M residual = %.6e deg\n', rad2deg(pmp.finalMError));
    fprintf('  status           = %s\n', pmp.status);

    lyap = simulate_lyapunov_slot(sat, slot, pmp, cfg);
    lyapResults(k).satId = sat.id;
    lyapResults(k).slotId = slot.id;
    lyapResults(k).solution = lyap;

    fprintf('Lyapunov solution:\n');
    fprintf('  final time       = %.6f days\n', lyap.tf/86400);
    fprintf('  arrival epoch    = %s\n', char(lyap.arrivalEpoch));
    fprintf('  delta-V          = %.6f m/s\n', lyap.deltaV_mps);
    fprintf('  final a error    = %.6e km\n', lyap.finalAError_km);
    fprintf('  final RAAN error = %.6e deg\n', rad2deg(lyap.finalRAANError));
    fprintf('  final M residual = %.6e deg\n', rad2deg(lyap.finalMError));
    fprintf('  status           = %s\n', lyap.status);
end

%% Build and save summary table
Sat = strings(numel(assignment),1);
Slot = strings(numel(assignment),1);
PMP_arrival_epoch = strings(numel(assignment),1);
LYAP_arrival_epoch = strings(numel(assignment),1);
PMP_tf_days = zeros(numel(assignment),1);
PMP_coast_days = zeros(numel(assignment),1);
PMP_raise_days = zeros(numel(assignment),1);
PMP_DV_mps = zeros(numel(assignment),1);
PMP_RAAN_err_deg = zeros(numel(assignment),1);
PMP_M_err_deg = zeros(numel(assignment),1);
LYAP_tf_days = zeros(numel(assignment),1);
LYAP_DV_mps = zeros(numel(assignment),1);
LYAP_RAAN_err_deg = zeros(numel(assignment),1);
LYAP_M_err_deg = zeros(numel(assignment),1);

for k = 1:numel(assignment)
    Sat(k) = string(pmpResults(k).satId);
    Slot(k) = string(pmpResults(k).slotId);

    pmp = pmpResults(k).solution;
    lyap = lyapResults(k).solution;

    PMP_arrival_epoch(k) = string(char(pmp.arrivalEpoch));
    LYAP_arrival_epoch(k) = string(char(lyap.arrivalEpoch));

    PMP_tf_days(k) = pmp.tf/86400;
    PMP_coast_days(k) = pmp.tCoast/86400;
    PMP_raise_days(k) = pmp.tRaise/86400;
    PMP_DV_mps(k) = pmp.deltaV_mps;
    PMP_RAAN_err_deg(k) = rad2deg(pmp.finalRAANError);
    PMP_M_err_deg(k) = rad2deg(pmp.finalMError);

    LYAP_tf_days(k) = lyap.tf/86400;
    LYAP_DV_mps(k) = lyap.deltaV_mps;
    LYAP_RAAN_err_deg(k) = rad2deg(lyap.finalRAANError);
    LYAP_M_err_deg(k) = rad2deg(lyap.finalMError);
end

summaryTable = table(Sat, Slot, PMP_arrival_epoch, PMP_tf_days, PMP_coast_days, ...
    PMP_raise_days, PMP_DV_mps, PMP_RAAN_err_deg, PMP_M_err_deg, ...
    LYAP_arrival_epoch, LYAP_tf_days, LYAP_DV_mps, LYAP_RAAN_err_deg, ...
    LYAP_M_err_deg);

resultsDir = fullfile(projectRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

save(fullfile(resultsDir, 'starlink_assignment_results.mat'), ...
    'cfg', 'assignment', 'C', 'pairEst', 'pmpResults', 'lyapResults', 'summaryTable');

writetable(summaryTable, fullfile(resultsDir, 'assignment_summary.csv'));

fprintf('\n================ SUMMARY TABLE ================\n');
disp(summaryTable);

fprintf('\nSaved results to:\n');
fprintf('  %s\n', fullfile(resultsDir, 'starlink_assignment_results.mat'));
fprintf('  %s\n', fullfile(resultsDir, 'assignment_summary.csv'));

%% Plots
plot_starlink_map(cfg, assignment);
plot_control_results(cfg, pmpResults, lyapResults);
