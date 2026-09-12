%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/06/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function figures = plot_interaction_chance_null_grid( ...
    nullData, observedMeans, fontSize)

assert(isequal(size(nullData),[3 4]), ...
    'nullData must be a 3-by-4 cell array.');

assert(isequal(size(observedMeans),[3 4]), ...
    'observedMeans must be a 3-by-4 numeric matrix.');

columnTitles = { ...
    'A-to-C Pre-Switch', ...
    'C-to-A Pre-Switch', ...
    'A-to-C Post-Switch', ...
    'C-to-A Post-Switch'};

metricLabels = { ...
    'Interaction frequency (Hz)', ...
    'Mean interaction duration (s)', ...
    'Total interaction time (s)'};

a2cColor = [1.0 0.4 0.4];
c2aColor = [0.2 0.6 1.0];

columnColors = { ...
    a2cColor, c2aColor, ...
    a2cColor, c2aColor};

figures = gobjects(12,1);
axesHandles = gobjects(3,4);

panelIndex = 0;

for metricIndex = 1:3

    % Same histogram bins within each metric row
    pooledNull = vertcat(nullData{metricIndex,:});
    pooledNull = pooledNull(isfinite(pooledNull));

    nullMinimum = min(pooledNull);
    nullMaximum = max(pooledNull);

    if nullMaximum <= nullMinimum
        nullMinimum = nullMinimum - 0.5;
        nullMaximum = nullMaximum + 0.5;
    end

    binEdges = linspace(nullMinimum,nullMaximum,51);

    % Same x-axis range within each metric row
    rowValues = [pooledNull; observedMeans(metricIndex,:)'];
    rowValues = rowValues(isfinite(rowValues));

    xMinimum = min(rowValues);
    xMaximum = max(rowValues);
    xPadding = 0.05*(xMaximum-xMinimum);

    if ~isfinite(xPadding) || xPadding <= 0
        xPadding = 0.5;
    end

    for columnIndex = 1:4

        panelIndex = panelIndex + 1;

        fig = figure( ...
            'Color','w', ...
            'Position',[680 530 420 330]);

        figures(panelIndex) = fig;

        ax = axes(fig);
        axesHandles(metricIndex,columnIndex) = ax;
        hold(ax,'on');

        nullValues = nullData{metricIndex,columnIndex};
        nullValues = nullValues(isfinite(nullValues));

        plotColor = columnColors{columnIndex};

        histogram(ax,nullValues,binEdges, ...
            'Normalization','probability', ...
            'FaceColor',plotColor, ...
            'FaceAlpha',0.45, ...
            'EdgeColor','none');

        chanceMean = mean(nullValues,'omitnan');
        observedMean = observedMeans(metricIndex,columnIndex);

        % Dashed line = shuffled chance mean
        xline(ax,chanceMean,'--', ...
            'Color',plotColor, ...
            'LineWidth',2);

        % Solid line = observed mean
        xline(ax,observedMean,'-', ...
            'Color',plotColor, ...
            'LineWidth',2.5);

        xlim(ax,[xMinimum-xPadding, xMaximum+xPadding]);

        title(ax,columnTitles{columnIndex}, ...
            'FontWeight','normal');

        xlabel(ax,metricLabels{metricIndex});
        ylabel(ax,'Probability');

        panelLetter = char('A' + panelIndex - 1);

        text(ax,-0.14,1.07,panelLetter, ...
            'Units','normalized', ...
            'FontWeight','bold', ...
            'FontSize',fontSize, ...
            'HorizontalAlignment','left');

        set(ax, ...
            'FontSize',max(fontSize-3,10), ...
            'LineWidth',0.8, ...
            'TickDir','in');

        box(ax,'off');
    end

    drawnow;

    rowYMaximum = max(arrayfun( ...
        @(currentAxes) currentAxes.YLim(2), ...
        axesHandles(metricIndex,:)));

    set(axesHandles(metricIndex,:), ...
        'YLim',[0 rowYMaximum]);

end

end
