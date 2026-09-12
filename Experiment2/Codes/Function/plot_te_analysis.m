%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_te_analysis(TE, visibleState, figureDir)

preColor = [0.72 0.62 0.85];
postColor = [0.36 0.18 0.55];
axisColor = [33 33 33] / 255;

requiredFields = {'PoolNames','NetTEByLocationPhase', ...
    'ChanceByLocationPhase','NullByLocationPhase', ...
    'TAU','DelayValues','DelaySweep'};
assert(all(isfield(TE,requiredFields)), ...
    'Results.TE does not contain all fields required for TE plotting.');

% Common-TAU selection
fig = make_figure(visibleState,[680 530 327*1.5 348]);
ax = axes(fig);

plot(ax,TE.DelayValues,TE.DelaySweep, ...
    '-o','Color',postColor,'MarkerFaceColor',preColor, ...
    'MarkerEdgeColor',postColor,'LineWidth',1.4,'MarkerSize',4);
hold(ax,'on');
xline(ax,TE.TAU,'--','Color',[0.45 0.45 0.45],'LineWidth',1.1);

xlabel(ax,'Delay, \tau (samples)','Interpreter','tex','FontSize',11);
ylabel(ax,'Mean robot-to-fish TE (bits)','Interpreter','tex','FontSize',11);
xticks(ax,TE.DelayValues);
set(ax,'FontSize',9,'LineWidth',0.5,'TickDir','in', ...
    'XColor',axisColor,'YColor',axisColor);
box(ax,'off');

save_figure(fig,figureDir,'TE_Delay_Selection');

% Shuffled-null distributions
nPools = size(TE.NullByLocationPhase,1);
assert(nPools == 2 && size(TE.NullByLocationPhase,2) == 2, ...
    'Expected TE null data organized as 2 Locations x 2 Phases.');

figurePosition = [680 530 327*1.5 348];
leftMargin = 0.14;
rightMargin = 0.04;
bottom = 0.23;
height = 0.70;
panelWidth = 1-leftMargin-rightMargin;
segmentGap = 0.035;
leftFraction = 0.50;
leftWidth = panelWidth * leftFraction;
rightWidth = panelWidth - leftWidth - segmentGap;
breakPoint = 0.01;

for g = 1:nPools
    fig = make_figure(visibleState,figurePosition);
    panelX = leftMargin;
    leftPosition = [panelX bottom leftWidth height];
    rightPosition = [panelX+leftWidth+segmentGap bottom rightWidth height];

    axLeft = axes(fig,'Position',leftPosition);
    axRight = axes(fig,'Position',rightPosition);
    hold(axLeft,'on');
    hold(axRight,'on');

    nullPre = TE.NullByLocationPhase{g,1}(:);
    nullPost = TE.NullByLocationPhase{g,2}(:);
    nullPre = nullPre(isfinite(nullPre));
    nullPost = nullPost(isfinite(nullPost));

    allNull = [nullPre;nullPost];
    assert(~isempty(allNull), ...
        'TE null data contain no finite values.');

    chancePre = TE.ChanceByLocationPhase{g,1};
    chancePost = TE.ChanceByLocationPhase{g,2};
    observedPre = mean(TE.NetTEByLocationPhase{g,1},'omitnan');
    observedPost = mean(TE.NetTEByLocationPhase{g,2},'omitnan');

    nullMin = min(allNull);
    nullMax = max(allNull);
    if nullMin == nullMax
        nullMin = nullMin - eps(max(abs(nullMin),1));
        nullMax = nullMax + eps(max(abs(nullMax),1));
    end
    nullEdges = linspace(nullMin,nullMax,41);

    histogram(axLeft,nullPre,nullEdges, ...
        'FaceColor',preColor,'FaceAlpha',0.55, ...
        'EdgeColor','none','Normalization','probability');
    histogram(axLeft,nullPost,nullEdges, ...
        'FaceColor',postColor,'FaceAlpha',0.55, ...
        'EdgeColor','none','Normalization','probability');

    xline(axLeft,chancePre,'--','Color',preColor,'LineWidth',1.4);
    xline(axLeft,chancePost,'--','Color',postColor,'LineWidth',1.4);
    xline(axRight,observedPre,'-','Color',preColor,'LineWidth',2.0);
    xline(axRight,observedPost,'-','Color',postColor,'LineWidth',2.0);

    xUpper = 0.01 * ceil( ...
        max([allNull;observedPre;observedPost]) / 0.01);
    xUpper = max(xUpper,breakPoint+0.01);

    preCounts = histcounts(nullPre,nullEdges, ...
        'Normalization','probability');
    postCounts = histcounts(nullPost,nullEdges, ...
        'Normalization','probability');
    yUpper = 1.08 * max([preCounts postCounts]);
    yUpper = max(yUpper,eps);

    leftData = [allNull;chancePre;chancePost];
    leftDataSpan = max(leftData)-min(leftData);
    leftPadding = max(0.05*leftDataSpan,0.0001);
    leftLower = max(0,min(leftData)-leftPadding);
    leftUpper = min(breakPoint,max(leftData)+leftPadding);

    xlim(axLeft,[leftLower leftUpper]);
    firstLeftTick = 0.001 * ceil(leftLower/0.001);
    lastLeftTick = 0.001 * floor(leftUpper/0.001);
    leftTicks = firstLeftTick:0.001:lastLeftTick;
    if isempty(leftTicks)
        leftTicks = [leftLower leftUpper];
    end
    xticks(axLeft,leftTicks);
    xticklabels(axLeft,compose('%.3f',leftTicks));
    ylim(axLeft,[0 yUpper]);

    xlim(axRight,[breakPoint xUpper]);
    xticks(axRight,(breakPoint+0.01):0.01:xUpper);
    ylim(axRight,[0 yUpper]);

    ylabel(axLeft,'Probability','Interpreter','tex','FontSize',10);

    text(axLeft,0.04,0.94,string(TE.PoolNames{g}), ...
        'Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontSize',10,'FontWeight','bold');

    set(axLeft,'FontSize',8,'LineWidth',0.5,'TickDir','in', ...
        'XColor',axisColor,'YColor',axisColor);
    set(axRight,'FontSize',8,'LineWidth',0.5,'TickDir','in', ...
        'XColor',axisColor,'YColor','none');
    xtickangle(axLeft,45);
    xtickangle(axRight,45);
    box(axLeft,'off');
    box(axRight,'off');

    annotation(fig,'textbox', ...
        [panelX bottom-0.20 panelWidth 0.06], ...
        'String','Robot-to-fish TE (bits)', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'EdgeColor','none','FontSize',10,'Interpreter','tex');

    slashWidth = 0.005;
    slashHeight = 0.012;
    leftBreakX = leftPosition(1) + leftPosition(3);
    rightBreakX = rightPosition(1);
    annotation(fig,'line', ...
        [leftBreakX-slashWidth leftBreakX+slashWidth], ...
        [bottom-slashHeight bottom+slashHeight], ...
        'Color',axisColor,'LineWidth',0.8);
    annotation(fig,'line', ...
        [rightBreakX-slashWidth rightBreakX+slashWidth], ...
        [bottom-slashHeight bottom+slashHeight], ...
        'Color',axisColor,'LineWidth',0.8);

    poolFileName = regexprep(char(string(TE.PoolNames{g})), ...
        '[^A-Za-z0-9_-]','_');
    save_figure(fig,figureDir, ...
        ['TE_Shuffled_Null_Distribution_' poolFileName]);
end

end
