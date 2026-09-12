%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_neuron_raster_axis(ax, timeMin, M, titleText, climValues)
imagesc(ax, timeMin, 1:size(M, 1), M);

set(ax, ...
    'YDir', 'normal', ...
    'FontSize', 11, ...
    'TickLabelInterpreter', 'none', ...
    'XColor', [33 33 33]/255, ...
    'YColor', [33 33 33]/255, ...
    'TickDir', 'in', ...
    'LineWidth', 0.5);
title(ax, titleText, 'Interpreter', 'none', 'FontWeight', 'normal');
xlim(ax, [0 10]);
ylim(ax, [0.5, max(size(M, 1), 1) + 0.5]);
caxis(ax, climValues);
box(ax, 'on');

end
