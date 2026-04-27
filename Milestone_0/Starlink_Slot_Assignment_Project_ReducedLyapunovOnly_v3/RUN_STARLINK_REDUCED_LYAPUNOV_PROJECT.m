clc; clear; close all;

%% ============================================================
% STARLINK-LIKE SLOT ASSIGNMENT PROJECT
% REDUCED LYAPUNOV ONLY - NO PMP (V3 PLOT FIX)
%
% This version intentionally removes all PMP logic.  It reproduces the
% Project-3 style plots using only the reduced Lyapunov dynamics:
%
%   x = [a; Delta_RAAN; Delta_M]
%
%   a_dot          = 2 u_T / n(a)
%   Delta_RAANdot = OmegaDot_J2(a) - OmegaDot_J2(aF)
%   Delta_Mdot    = n(a) - n(aF)
%
% Controller:
%
%   V = 1/2 (a-aF)^2
%   u_T = sat[-kA(a-aF), +/- uMax]
%
% RAAN targeting is handled by a coast phase at the lower parking/drift
% altitude before the Lyapunov controlled raise.
%% ============================================================

thisFile = mfilename('fullpath');
projectRoot = fileparts(thisFile);
addpath(genpath(projectRoot));

cfg = create_starlink_reduced_lyapunov_scenario();

fprintf('\n============================================================\n');
fprintf(' STARLINK-LIKE SLOT ASSIGNMENT - REDUCED LYAPUNOV ONLY\n');
fprintf('============================================================\n');
fprintf('Launch epoch                  : %s\n', char(cfg.launchEpoch));
fprintf('Parking altitude              : %.3f km\n', cfg.hParking_km);
fprintf('Final shell altitude          : %.3f km\n', cfg.hFinal_km);
fprintf('Inclination                   : %.4f deg\n', cfg.inc_deg);
fprintf('Tangential acceleration limit : %.3e km/s^2 = %.4f mm/s^2\n', ...
    cfg.uMax_km_s2, 1e6*cfg.uMax_km_s2);
fprintf('Number of satellites          : %d\n', numel(cfg.sats));
fprintf('Number of slots               : %d\n', numel(cfg.slots));
fprintf('Controller                    : reduced Lyapunov only, no PMP\n');
fprintf('Delta-V assignment margin     : %.1f %%\n', 100*cfg.assignment.dvMarginFraction);

rateP = j2_raan_rate(cfg.aParking, cfg.inc_rad, cfg.Earth);
rateF = j2_raan_rate(cfg.aFinal,   cfg.inc_rad, cfg.Earth);
fprintf('\nJ2 RAAN drift at parking      : %.6f deg/day\n', rad2deg(rateP)*86400);
fprintf('J2 RAAN drift at final shell  : %.6f deg/day\n', rad2deg(rateF)*86400);
fprintf('Differential RAAN drift       : %.6f deg/day\n', rad2deg(rateP-rateF)*86400);

[Tr, dOmR, dMR] = compute_raise_integrals(cfg.aParking, cfg.aFinal, cfg.uMax_km_s2, cfg);
fprintf('\nNominal saturated Lyapunov raise time : %.6f days\n', Tr/86400);
fprintf('RAAN gained during raise relative to final shell : %.6f deg\n', rad2deg(dOmR));
fprintf('Mean anomaly gained during raise relative to final shell : %.6f deg\n', mod(rad2deg(dMR),360));

%% Initial satellites
fprintf('\n================ INITIAL SATELLITES ================\n');
for i = 1:numel(cfg.sats)
    s = cfg.sats(i);
    fprintf('%s: a0 = %.3f km, altitude = %.3f km, RAAN0 = %.6f deg, M0 = %.6f deg\n', ...
        s.id, s.a0, s.a0 - cfg.Earth.radius, s.raan0_deg, s.M0_deg);
end

%% Open slots
fprintf('\n================ OPEN SLOTS ================\n');
for j = 1:numel(cfg.slots)
    q = cfg.slots(j);
    fprintf('%s: aF = %.3f km, altitude = %.3f km, RAAN0 = %.6f deg, M0 = %.6f deg, nominal epoch = %s\n', ...
        q.id, q.aF, q.aF - cfg.Earth.radius, q.raan0_deg, q.M0_deg, char(q.nominalEpoch));
end

%% Evaluate every satellite-slot pair
fprintf('\n================ EVALUATING ALL REDUCED-LYAPUNOV PAIRS ================\n');
[DV, solGrid, estGrid] = evaluate_reduced_lyapunov_delta_v_matrix(cfg);

%% Print matrix
print_reduced_lyapunov_delta_v_matrix(DV, cfg);

%% Assignment search
[best, allCases] = find_best_assignment_reduced_lyapunov(DV, cfg);
print_best_reduced_lyapunov_assignment(best, cfg, allCases);

%% Detailed report
fprintf('\n================ BEST SCENARIO DETAILED TRANSFER REPORT ================\n');
for k = 1:numel(best.assignment)
    aBest = best.assignment(k);
    sat = cfg.sats(aBest.satIdx);
    slot = cfg.slots(aBest.slotIdx);
    sol = solGrid(aBest.satIdx, aBest.slotIdx);

    fprintf('\n--- Best pair %d: %s -> %s ---\n', k, sat.id, slot.id);
    fprintf('  coast time       = %9.4f days\n', sol.tCoast/86400);
    fprintf('  controlled time  = %9.4f days\n', sol.tControlled/86400);
    fprintf('  final time       = %9.4f days\n', sol.tf/86400);
    fprintf('  arrival epoch    = %s\n', char(sol.arrivalEpoch));
    fprintf('  max |u_T|        = %9.4f mm/s^2\n', 1e6*max(abs(sol.u)));
    fprintf('  raw Lyapunov DV  = %9.4f m/s\n', sol.rawDeltaV_mps);
    fprintf('  cleanup DV       = %9.4f m/s\n', sol.cleanupDeltaV_mps);
    fprintf('  total DV         = %9.4f m/s\n', sol.deltaV_mps);
    fprintf('  final a error    = %+9.4e km\n', sol.finalAError_km);
    fprintf('  final RAAN error = %+9.5f deg\n', rad2deg(sol.finalRAANError));
    fprintf('  final M residual = %+9.5f deg\n', rad2deg(sol.finalMError));
    fprintf('  minimum altitude = %9.3f km\n', sol.minAltitude_km);
    fprintf('  status           = %s\n', char(sol.status));
end

%% Save outputs
resultsDir = fullfile(projectRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

summaryTable = build_reduced_lyapunov_summary_table(cfg, best, solGrid);

save(fullfile(resultsDir, 'reduced_lyapunov_only_results.mat'), ...
    'cfg', 'DV', 'solGrid', 'estGrid', 'best', 'allCases', 'summaryTable');

writetable(allCases, fullfile(resultsDir, 'reduced_lyapunov_assignment_permutations.csv'));
writetable(summaryTable, fullfile(resultsDir, 'reduced_lyapunov_best_assignment_summary.csv'));
writematrix_with_labels(DV.total_mps, cfg, fullfile(resultsDir, 'delta_v_matrix_reduced_lyapunov.csv'));

fprintf('\n================ BEST ASSIGNMENT SUMMARY TABLE ================\n');
disp(summaryTable);

fprintf('\nSaved results to:\n');
fprintf('  %s\n', fullfile(resultsDir, 'reduced_lyapunov_only_results.mat'));
fprintf('  %s\n', fullfile(resultsDir, 'reduced_lyapunov_assignment_permutations.csv'));
fprintf('  %s\n', fullfile(resultsDir, 'reduced_lyapunov_best_assignment_summary.csv'));
fprintf('  %s\n', fullfile(resultsDir, 'delta_v_matrix_reduced_lyapunov.csv'));

%% Plots: Project-3 style plus the requested summary plots
plot_starlink_map_reduced_lyapunov(cfg, best.assignment);
plot_delta_v_matrix_reduced_lyapunov(DV, cfg);
plot_reduced_lyapunov_histories(cfg, best, solGrid);
plot_best_orbits_3d_reduced_lyapunov(cfg, best, solGrid);
plot_sma_safety_check_reduced_lyapunov(cfg, best, solGrid);
plot_control_timeline_reduced_lyapunov(cfg, best, solGrid);

save_all_open_figures(resultsDir);

fprintf('\nAll plots generated and saved to:\n  %s\n', resultsDir);
