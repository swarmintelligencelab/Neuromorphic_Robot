%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_timebin_trajectories(metricCells, conditionDefs, colors, yLabel, plotTitle, ...
    visibleState, figureDir, fileBase)

nConditions = numel(metricCells);
nBins = size(metricCells{1}, 1);
time = 1:nBins;

fig = make_figure(visibleState, [680 530 327*1.5 348]);
set(fig, 'Position', [680 530 327*1.5 348]);
hold on;

for c = 1:nConditions
    M = metricCells{c};
    mu = mean(M, 2, 'omitnan');
    sem = sem_columns(M);
    errorbar(time, mu, sem, 'o-', 'Color', colors(c, :), ...
        'LineWidth', 2, 'MarkerSize', 5, 'CapSize', 4, ...
        'MarkerFaceColor', colors(c, :));
end

xlim([1, nBins]);
xlabel('Time bin');
ylabel(yLabel);
title(plotTitle);
legend({conditionDefs.display}, 'Location', 'bestoutside', 'Interpreter', 'none');
set(gca, ...
    'FontSize', 16, ...
    'TickLabelInterpreter', 'none', ...
    'XColor', [33 33 33]/255, ...
    'YColor', [33 33 33]/255, ...
    'TickDir', 'in', ...
    'LineWidth', 0.5);
set(gcf, 'Color', 'w');
grid off;
hold off;

save_figure(fig, figureDir, fileBase);

end
