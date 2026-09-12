%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_neuron_raster_pair(leftByCondition, binaryByCondition, conditionDefs, plotTimeMin, ...
    leftCLim, leftTitle, leftColorbarLabel, superTitle, visibleState, figureDir, fileBase)

nConditions = numel(leftByCondition);
fig = make_figure(visibleState, [100 100 1250 1050]);
tiledlayout(nConditions, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for c = 1:nConditions
    axLeft = nexttile((c - 1) * 2 + 1);
    plot_neuron_raster_axis(axLeft, plotTimeMin, leftByCondition{c}, ...
        conditionDefs(c).display, leftCLim);
    colormap(axLeft, neuron_raw_colormap());
    if c == 1
        title(axLeft, leftTitle, 'FontWeight', 'normal');
    end
    if c == nConditions
        xlabel(axLeft, 'Time aligned to switch (min)');
    end
    ylabel(axLeft, 'Trial index');
    xline(axLeft, 5, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);

    cbLeft = colorbar(axLeft);
    cbLeft.Label.String = leftColorbarLabel;
    cbLeft.Label.Interpreter = 'tex';

    axBin = nexttile((c - 1) * 2 + 2);
    plot_neuron_raster_axis(axBin, plotTimeMin, binaryByCondition{c}, ...
        conditionDefs(c).display, [0 1]);
    colormap(axBin, neuron_binary_colormap());
    if c == 1
        title(axBin, 'Binary firing', 'FontWeight', 'normal');
    end
    if c == nConditions
        xlabel(axBin, 'Time aligned to switch (min)');
    end
    ylabel(axBin, 'Trial index');
    xline(axBin, 5, '--k', 'LineWidth', 1.0);

    cbBin = colorbar(axBin);
    cbBin.Ticks = [0 1];
end

sgtitle(superTitle, 'Interpreter', 'none');
save_figure(fig, figureDir, fileBase);

end
