function plot_delta_v_matrix_reduced_lyapunov(DV, cfg)
%PLOT_DELTA_V_MATRIX_REDUCED_LYAPUNOV Heatmap using imagesc.

    figure('Color','w','Name','plot_delta_v_matrix_reduced_lyapunov');
    imagesc(DV.total_mps);
    colorbar;
    axis equal tight;
    set(gca, 'XTick', 1:numel(cfg.slots), 'XTickLabel', {cfg.slots.id});
    set(gca, 'YTick', 1:numel(cfg.sats),  'YTickLabel', {cfg.sats.id});
    xlabel('Slots');
    ylabel('Satellites');
    title('Reduced Lyapunov Total \DeltaV Matrix [m/s]');

    for i = 1:numel(cfg.sats)
        for j = 1:numel(cfg.slots)
            text(j, i, sprintf('%.2f', DV.total_mps(i,j)), ...
                'HorizontalAlignment','center', 'Color','w', 'FontWeight','bold');
        end
    end
end
