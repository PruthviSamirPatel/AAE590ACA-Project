function plot_control_timeline_reduced_lyapunov(cfg, best, solGrid)
%PLOT_CONTROL_TIMELINE_REDUCED_LYAPUNOV Coast/control timeline.
%
% This plot intentionally shows a hard two-phase switch:
%   mode = 0 : coast/hold at parking orbit for J2 RAAN drift
%   mode = 1 : reduced Lyapunov tangential thrusting to raise SMA
%
% The vertical dotted line is the analytically planned coast time.

    figure('Color','w','Name','control_timeline');
    n = numel(best.assignment);

    for k = 1:n
        i = best.assignment(k).satIdx;
        j = best.assignment(k).slotIdx;
        sol = solGrid(i,j);

        subplot(n,1,k);
        stairs(sol.t/86400, sol.mode(:), 'LineWidth',1.6);
        hold on; grid on; box on;
        xline(sol.tCoast/86400, ':', 'LineWidth',1.1);

        ylim([-0.15 1.15]);
        yticks([0 1]);
        yticklabels({'coast','thrust'});
        ylabel(sprintf('%s->%s', cfg.sats(i).id, cfg.slots(j).id));

        if k == 1
            title('Reduced Lyapunov Timeline: coast for J2 drift, then thrust');
        end
    end

    xlabel('Time since launch [days]');
end
