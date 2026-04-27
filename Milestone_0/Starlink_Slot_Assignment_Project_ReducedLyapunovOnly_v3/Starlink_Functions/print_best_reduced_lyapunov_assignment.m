function print_best_reduced_lyapunov_assignment(best, cfg, allCases)
%PRINT_BEST_REDUCED_LYAPUNOV_ASSIGNMENT Console summary.

    fprintf('\n================ ECONOMICAL REDUCED-LYAPUNOV ASSIGNMENT ================\n');
    fprintf('Objective for each selected transfer: (1 + %.2f)*DeltaV = %.2f*DeltaV\n', ...
        cfg.assignment.dvMarginFraction, 1 + cfg.assignment.dvMarginFraction);

    fprintf('\n%8s %8s %16s %16s\n', 'Sat', 'Slot', 'DV [m/s]', 'DV+20% [m/s]');
    for k = 1:numel(best.assignment)
        a = best.assignment(k);
        fprintf('%8s %8s %16.4f %16.4f\n', ...
            cfg.sats(a.satIdx).id, cfg.slots(a.slotIdx).id, ...
            a.deltaV_mps, a.deltaV_withMargin_mps);
    end

    fprintf('\nTotal raw Delta-V for best Lyapunov scenario : %.6f m/s\n', best.rawDeltaV_mps);
    fprintf('Total objective Delta-V with 20%% margin      : %.6f m/s\n', best.objectiveDeltaV_mps);

    fprintf('\n--- Top 10 Lyapunov assignments by objective ---\n');
    disp(allCases(1:min(10,height(allCases)),:));
end
