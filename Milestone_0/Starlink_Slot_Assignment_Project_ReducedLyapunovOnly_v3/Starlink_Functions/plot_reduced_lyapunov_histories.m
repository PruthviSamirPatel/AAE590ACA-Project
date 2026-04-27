function plot_reduced_lyapunov_histories(cfg, best, solGrid)
%PLOT_REDUCED_LYAPUNOV_HISTORIES Plot selected transfer histories.
%
% Notes:
%   - The control command is plotted in mm/s^2 so the Project-3 low-thrust
%     value of 0.2 mm/s^2 is visible.  Plotting in m/s^2 makes SAT01/SAT02
%     look like zero because 0.0002 m/s^2 is tiny on MATLAB's default axis.
%   - Coast-to-raise switching is intentionally discontinuous in this
%     reduced two-phase model: first wait for J2 RAAN drift, then thrust.

    n = numel(best.assignment);

    %% Altitude histories
    figure('Color','w','Name','altitude_histories');
    for k = 1:n
        subplot(n,1,k);
        sol = get_best_sol(k, best, solGrid);
        plot(sol.t/86400, sol.a - cfg.Earth.radius, 'LineWidth',1.5); hold on;
        yline(cfg.hFinal_km, 'k:');
        yline(cfg.hParking_km, 'm:');
        grid on;
        ylabel(sprintf('%s->%s\nh [km]', sol.satId, sol.slotId));
        if k == 1
            title('Reduced Lyapunov Altitude Histories');
            legend('Transfer','Final shell','Parking shell','Location','best');
        end
    end
    xlabel('Time since launch [days]');

    %% RAAN error histories
    figure('Color','w','Name','raan_error_histories');
    for k = 1:n
        subplot(n,1,k);
        sol = get_best_sol(k, best, solGrid);
        plot(sol.t/86400, rad2deg(sol.relRAAN), 'LineWidth',1.5); hold on;
        yline(0, 'k:');
        grid on;
        ylabel(sprintf('%s->%s\n\\Delta\\Omega [deg]', sol.satId, sol.slotId));
        if k == 1
            title('Reduced Lyapunov RAAN Error Histories');
        end
    end
    xlabel('Time since launch [days]');

    %% Mean-anomaly residual histories
    figure('Color','w','Name','mean_anomaly_histories');
    for k = 1:n
        subplot(n,1,k);
        sol = get_best_sol(k, best, solGrid);
        plot(sol.t/86400, rad2deg(sol.relM), 'LineWidth',1.5); hold on;
        yline(0, 'k:');
        grid on;
        ylabel(sprintf('%s->%s\n\\Delta M [deg]', sol.satId, sol.slotId));
        if k == 1
            title('Reduced Lyapunov Mean-Anomaly Residual Histories');
        end
    end
    xlabel('Time since launch [days]');

    %% Tangential control histories
    figure('Color','w','Name','control_histories');
    uMax_mm_s2 = 1e6*cfg.uMax_km_s2;
    for k = 1:n
        subplot(n,1,k);
        sol = get_best_sol(k, best, solGrid);

        stairs(sol.t/86400, 1e6*sol.u, 'LineWidth',1.4); hold on;
        yline( uMax_mm_s2, 'k:', 'LineWidth',0.8);
        yline(0, 'k-', 'LineWidth',0.5);
        grid on;
        ylabel(sprintf('%s->%s\nu_T [mm/s^2]', sol.satId, sol.slotId));
        ylim([-0.05*uMax_mm_s2, 1.20*uMax_mm_s2]);

        if k == 1
            title('Reduced Lyapunov Tangential Control Histories');
            legend('u_T','u_{max}','zero','Location','best');
        end
    end
    xlabel('Time since launch [days]');

    %% Cumulative raw Delta-V histories
    figure('Color','w','Name','cumulative_dv_histories');
    for k = 1:n
        subplot(n,1,k);
        sol = get_best_sol(k, best, solGrid);
        cumDV = cumtrapz(sol.t, abs(sol.u))*1000;
        plot(sol.t/86400, cumDV, 'LineWidth',1.5);
        grid on;
        ylabel(sprintf('%s->%s\nDV [m/s]', sol.satId, sol.slotId));
        if k == 1
            title('Cumulative Raw Lyapunov \DeltaV Histories');
        end
    end
    xlabel('Time since launch [days]');
end

function sol = get_best_sol(k, best, solGrid)
    i = best.assignment(k).satIdx;
    j = best.assignment(k).slotIdx;
    sol = solGrid(i,j);
end
