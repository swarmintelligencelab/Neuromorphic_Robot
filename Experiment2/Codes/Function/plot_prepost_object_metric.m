%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_prepost_object_metric(metricCells, yLabel, ~, ...
    visibleState, figureDir, fileBase, varargin)

chanceStatsFile = '';
if ~isempty(varargin)
    chanceStatsFile = varargin{1};
end

idxIn = [3 2 1];  % Rod, Rectangle, Fish model
idxOut = [6 5 4]; % Rod, Rectangle, Fish model
idx = [idxIn idxOut];

xCenters = [1 2 3 5 6 7];
objectNames = {'Rod','Rectangle','Fish model', ...
               'Rod','Rectangle','Fish model'};

preColor = [0.72 0.62 0.85];  % light purple
postColor = [0.36 0.18 0.55]; % dark purple

Y = nan(6,2);
SEM = nan(6,2);

for j = 1:6
    D = metricCells{idx(j)};
    assert(size(D,1) == 2, ...
        'Each metric matrix must contain two rows: pre-switch and post-switch.');

    for phase = 1:2
        values = D(phase,:);
        nValid = sum(~isnan(values));
        Y(j,phase) = mean(values,'omitnan');
        if nValid > 0
            SEM(j,phase) = std(values,0,'omitnan') / sqrt(nValid);
        end
    end
end

fig = make_figure(visibleState, [680 530 327*1.5 348]);
set(fig, 'Color', 'w', 'Position', [680 530 327*1.5 348]);
ax = axes(fig);
hold(ax,'on');

b = bar(ax,xCenters,Y,'grouped','BarWidth',0.75);
b(1).FaceColor = preColor;
b(1).EdgeColor = 'none';
b(1).LineWidth = 0.8;
b(2).FaceColor = postColor;
b(2).EdgeColor = 'none';
b(2).LineWidth = 0.8;

for phase = 1:2
    errorbar(ax,b(phase).XEndPoints,Y(:,phase),SEM(:,phase), ...
        'k','LineStyle','none','LineWidth',1.1,'CapSize',5);
end

xline(ax,4,'--','Color',[0.65 0.65 0.65],'LineWidth',1.1);
xlim(ax,[0.3 7.7]);
xticks(ax,xCenters);
xticklabels(ax,objectNames);
xtickangle(ax,30);
ylabel(ax,yLabel,'Interpreter','tex','FontSize',11);

originalYTicks = yticks(ax);

% Add space below y = 0 for the location labels
yl = ylim(ax);
yr = diff(yl);
ylim(ax,[yl(1)-0.14*yr, yl(2)]);

yticks(ax,originalYTicks(originalYTicks >= 0));

% Position the labels
yl = ylim(ax);
yText = yl(1) + 0.025*diff(yl);

text(ax,2,yText,'Inside tank', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',10);

text(ax,6,yText,'Outside tank', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',10);

if ~isempty(chanceStatsFile) && isfile(chanceStatsFile)
    nChanceMarkers = add_te_chance_significance_markers2( ...
        ax, chanceStatsFile, b, Y, SEM);
    fprintf('%s: added %d TE chance comparison marker(s) from %s\n', ...
        fileBase, nChanceMarkers, chanceStatsFile);
end

oldAxisColor = [33 33 33] / 255;

set(ax, ...
    'FontSize', 9, ...
    'LineWidth', 0.5, ...
    'TickDir', 'in', ...
    'XColor', oldAxisColor, ...
    'YColor', oldAxisColor);
box(ax,'off');

save_figure(fig, figureDir, fileBase);

end
