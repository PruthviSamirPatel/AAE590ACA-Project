function plot_control_results(cfg, pmpResults, lyapResults)
%PLOT_CONTROL_RESULTS Plot PMP and Lyapunov histories.

    n = numel(pmpResults);

    figure('Color','w');
    tiledlayout(n,1,'TileSpacing','compact');
    for k = 1:n
        nexttile;
        p = pmpResults(k).solution;
        l = lyapResults(k).solution;
        plot(p.t/86400, p.a - cfg.Earth.radius, 'LineWidth',1.5); hold on;
        plot(l.t/86400, l.a - cfg.Earth.radius, '--', 'LineWidth',1.5);
        yline(cfg.hFinal_km, 'k:');
        grid on;
        ylabel(sprintf('%s->%s\nh [km]', p.satId, p.slotId));
        if k == 1
            title('Altitude / SMA Histories');
            legend('PMP','Lyapunov','Final shell','Location','best');
        end
    end
    xlabel('Time since launch [days]');

    figure('Color','w');
    tiledlayout(n,1,'TileSpacing','compact');
    for k = 1:n
        nexttile;
        p = pmpResults(k).solution;
        l = lyapResults(k).solution;
        plot(p.t/86400, rad2deg(p.relRAAN), 'LineWidth',1.5); hold on;
        plot(l.t/86400, rad2deg(l.relRAAN), '--', 'LineWidth',1.5);
        yline(0, 'k:');
        grid on;
        ylabel(sprintf('%s->%s\n\\Delta\\Omega [deg]', p.satId, p.slotId));
        if k == 1
            title('RAAN Error Histories');
            legend('PMP','Lyapunov','Target','Location','best');
        end
    end
    xlabel('Time since launch [days]');

    figure('Color','w');
    tiledlayout(n,1,'TileSpacing','compact');
    for k = 1:n
        nexttile;
        p = pmpResults(k).solution;
        l = lyapResults(k).solution;
        plot(p.t/86400, rad2deg(p.relM), 'LineWidth',1.5); hold on;
        plot(l.t/86400, rad2deg(l.relM), '--', 'LineWidth',1.5);
        yline(0, 'k:');
        grid on;
        ylabel(sprintf('%s->%s\n\\Delta M [deg]', p.satId, p.slotId));
        if k == 1
            title('Mean-Anomaly Residual Histories');
            legend('PMP','Lyapunov','Target','Location','best');
        end
    end
    xlabel('Time since launch [days]');

    figure('Color','w');
    tiledlayout(n,1,'TileSpacing','compact');
    for k = 1:n
        nexttile;
        p = pmpResults(k).solution;
        l = lyapResults(k).solution;
        stairs(p.t/86400, 1e6*p.u, 'LineWidth',1.3); hold on;
        plot(l.t/86400, 1e6*l.u, '--', 'LineWidth',1.3);
        grid on;
        ylabel(sprintf('%s->%s\nu_T [mm/s^2]', p.satId, p.slotId));
        if k == 1
            title('Tangential Control Histories');
            legend('PMP','Lyapunov','Location','best');
        end
    end
    xlabel('Time since launch [days]');
end
