%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function nMarkers = add_te_chance_significance_markers( ...
    ax, statsFile, barHandles, meanValues, semValues)

T = readtable(statsFile, ...
    'FileType','text', ...
    'Delimiter',',', ...
    'ReadVariableNames',true, ...
    'VariableNamingRule','preserve');

requiredColumns = {'RobotShape','Location','Time','p_FDR'};
if ~all(ismember(requiredColumns,T.Properties.VariableNames))
    warning('TE chance statistics file has incompatible columns: %s',statsFile);
    nMarkers = 0;
    return;
end

yl = ylim(ax);
yRange = diff(yl);
if ~isfinite(yRange) || yRange <= 0
    yRange = 1;
end

nMarkers = 0;
maxMarkerY = yl(2);

for row = 1:height(T)
    pFdr = double(T.p_FDR(row));
    if ~isfinite(pFdr) || pFdr >= 0.05
        continue;
    end

    location = lower(strtrim(string(T.Location(row))));
    shape = lower(strtrim(string(T.RobotShape(row))));
    time = lower(strtrim(string(T.Time(row))));

    if location == "in"
        locationOffset = 0;
    elseif location == "out"
        locationOffset = 3;
    else
        warning('Could not map TE chance location: %s',location);
        continue;
    end

    if shape == "rod"
        conditionIndex = locationOffset + 1;
    elseif shape == "rectangle"
        conditionIndex = locationOffset + 2;
    elseif shape == "fish" || shape == "fish model" || shape == "fishmodel"
        conditionIndex = locationOffset + 3;
    else
        warning('Could not map TE chance stimulus morphology: %s',shape);
        continue;
    end

    if time == "pre" || contains(time,"pre-switch")
        phaseIndex = 1;
    elseif time == "post" || contains(time,"post-switch")
        phaseIndex = 2;
    else
        warning('Could not map TE chance phase: %s',time);
        continue;
    end

    markerX = barHandles(phaseIndex).XEndPoints(conditionIndex);
    markerY = meanValues(conditionIndex,phaseIndex) + ...
        semValues(conditionIndex,phaseIndex) + 0.025*yRange;

    text(ax,markerX,markerY,'#', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',13, ...
        'FontWeight','bold', ...
        'Interpreter','none', ...
        'HandleVisibility','off');

    maxMarkerY = max(maxMarkerY,markerY);
    nMarkers = nMarkers + 1;
end

if nMarkers > 0
    ylim(ax,[yl(1),max(yl(2),maxMarkerY + 0.05*yRange)]);
end

end
