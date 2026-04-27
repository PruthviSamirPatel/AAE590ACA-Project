function plot_starlink_map(cfg, assignment)
%PLOT_STARLINK_MAP Plot launch satellites and open slots in RAAN/M map.

    figure('Color','w');
    hold on; grid on; box on;

    satRAAN = [cfg.sats.raan0_deg];
    satM = [cfg.sats.M0_deg];

    slotRAAN = [cfg.slots.raan0_deg];
    slotM = [cfg.slots.M0_deg];

    scatter(satRAAN, satM, 80, 'filled', 'DisplayName','Launched satellites');
    scatter(slotRAAN, slotM, 100, 's', 'LineWidth',1.5, 'DisplayName','Open slots');

    for i = 1:numel(cfg.sats)
        text(satRAAN(i)+0.03, satM(i), cfg.sats(i).id);
    end

    for j = 1:numel(cfg.slots)
        text(slotRAAN(j)+0.03, slotM(j), cfg.slots(j).id);
    end

    for k = 1:numel(assignment)
        i = assignment(k).satIdx;
        j = assignment(k).slotIdx;
        plot([satRAAN(i), slotRAAN(j)], [satM(i), slotM(j)], 'k--', ...
             'HandleVisibility','off');
    end

    xlabel('RAAN at launch epoch [deg]');
    ylabel('Mean anomaly at launch epoch [deg]');
    title('RAAN / Mean-Anomaly Slot Assignment Map');
    legend('Location','best');
end
