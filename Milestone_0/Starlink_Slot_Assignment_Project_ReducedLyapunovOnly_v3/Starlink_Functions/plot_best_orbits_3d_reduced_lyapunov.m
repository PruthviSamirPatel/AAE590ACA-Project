function plot_best_orbits_3d_reduced_lyapunov(cfg, best, solGrid)
%PLOT_BEST_ORBITS_3D_REDUCED_LYAPUNOV 3-D Earth and selected transfer arcs.

    figure('Color','w','Name','plot_best_orbits_3d_reduced_lyapunov');
    hold on; grid on; axis equal;

    [xe, ye, ze] = sphere(60);
    surf(cfg.Earth.radius*xe, cfg.Earth.radius*ye, cfg.Earth.radius*ze, ...
        'FaceAlpha',0.45, 'EdgeColor','none', 'DisplayName','Earth');

    colors = lines(numel(best.assignment));

    for k = 1:numel(best.assignment)
        i = best.assignment(k).satIdx;
        j = best.assignment(k).slotIdx;
        sat = cfg.sats(i);
        slot = cfg.slots(j);
        sol = solGrid(i,j);

        R = trajectory_positions_from_solution(sol, sat, slot, cfg);

        plot3(R(1,:), R(2,:), R(3,:), 'LineWidth',1.5, 'Color',colors(k,:), ...
            'DisplayName',sprintf('%s to %s', sat.id, slot.id));

        % Initial parking orbit at launch epoch.
        th = linspace(0, 2*pi, 240);
        Ri = zeros(3,numel(th));
        for q = 1:numel(th)
            Ri(:,q) = circular_elements_to_rv(sat.a0, cfg.inc_rad, ...
                deg2rad(sat.raan0_deg), th(q));
        end
        plot3(Ri(1,:), Ri(2,:), Ri(3,:), ':', 'Color',colors(k,:), ...
            'HandleVisibility','off');

        % Final slot orbit at arrival epoch.
        OmDotF = j2_raan_rate(slot.aF, cfg.inc_rad, cfg.Earth);
        raanArrival = deg2rad(slot.raan0_deg) + OmDotF*sol.tf;
        Rf = zeros(3,numel(th));
        for q = 1:numel(th)
            Rf(:,q) = circular_elements_to_rv(slot.aF, cfg.inc_rad, raanArrival, th(q));
        end
        plot3(Rf(1,:), Rf(2,:), Rf(3,:), '--', 'Color',colors(k,:), ...
            'HandleVisibility','off');
    end

    xlabel('x [km]');
    ylabel('y [km]');
    zlabel('z [km]');
    title('Reduced Lyapunov Transfers: 3-D Earth, Parking Orbits, and Arrival Orbits');
    legend('Location','bestoutside');
    view(3);
end
