function plot_sma_safety_check_reduced_lyapunov(cfg, best, solGrid)
%PLOT_SMA_SAFETY_CHECK_REDUCED_LYAPUNOV
% Safety check in semi-major axis, not a separate altitude-only plot.
%
% The requirement h >= 200 km is equivalent to
%
%     a >= R_E + 200 km
%
% for this circular reduced model.  Separate subplots are used so SAT01 and
% SAT02 cannot hide each other when their curves overlap.

    hLimit = 200.0;                         % km
    aSafe  = cfg.Earth.radius + hLimit;     % km

    figure('Color','w','Name','sma_safety_check');
    n = numel(best.assignment);

    fprintf('\n================ SMA / 200 km SAFETY CHECK ================\n');

    for k = 1:n
        i = best.assignment(k).satIdx;
        j = best.assignment(k).slotIdx;
        sol = solGrid(i,j);

        subplot(n,1,k);
        plot(sol.t/86400, sol.a, 'LineWidth',1.5, ...
            'DisplayName',sprintf('%s to %s', sol.satId, sol.slotId));
        hold on; grid on; box on;

        yline(cfg.aFinal, 'k:', 'LineWidth',1.0, 'DisplayName','Final SMA');
        yline(cfg.aParking, 'm:', 'LineWidth',1.0, 'DisplayName','Parking SMA');
        yline(aSafe, 'r--', 'LineWidth',1.0, 'DisplayName','200 km safety SMA');

        ylabel(sprintf('%s->%s\na [km]', sol.satId, sol.slotId));

        if k == 1
            title('Reduced Lyapunov SMA Safety Check');
            legend('Location','best');
        end

        minA = min(sol.a);
        minH = minA - cfg.Earth.radius;
        if minA >= aSafe
            status = 'OK';
        else
            status = 'VIOLATION';
        end

        fprintf('%s -> %s: min a = %.3f km, min altitude = %.3f km, safety a = %.3f km    %s\n', ...
            cfg.sats(i).id, cfg.slots(j).id, minA, minH, aSafe, status);
    end

    xlabel('Time since launch [days]');
end
