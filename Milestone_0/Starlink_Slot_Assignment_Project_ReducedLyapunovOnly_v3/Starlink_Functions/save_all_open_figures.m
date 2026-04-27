function save_all_open_figures(resultsDir)
%SAVE_ALL_OPEN_FIGURES Save all open figures to PNG and FIG files.

    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        return
    end

    [~, order] = sort([figs.Number]);
    figs = figs(order);

    for k = 1:numel(figs)
        fig = figs(k);
        name = get(fig, 'Name');
        if isempty(name)
            name = sprintf('figure_%02d', fig.Number);
        end
        safeName = regexprep(char(name), '[^a-zA-Z0-9_\-]', '_');
        if isempty(safeName)
            safeName = sprintf('figure_%02d', fig.Number);
        end

        pngFile = fullfile(resultsDir, sprintf('%02d_%s.png', k, safeName));
        figFile = fullfile(resultsDir, sprintf('%02d_%s.fig', k, safeName));

        try
            saveas(fig, pngFile);
            savefig(fig, figFile);
        catch
            % Older MATLAB versions may not support savefig.
            saveas(fig, pngFile);
        end
    end
end
