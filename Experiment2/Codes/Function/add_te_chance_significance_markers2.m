%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function nMarkers = add_te_chance_significance_markers2( ...
    ax, statsFile, barHandles, meanValues, semValues)

% Read R chance-comparison results
T = readtable( ...
    statsFile, ...
    'FileType','text', ...
    'Delimiter',',', ...
    'ReadVariableNames',true, ...
    'VariableNamingRule','preserve');

requiredColumns = { ...
    'RobotShape', ...
    'Location', ...
    'Time', ...
    'ChanceMean', ...
    'p_FDR'};

if ~all(ismember(requiredColumns,T.Properties.VariableNames))

    warning( ...
        'Chance statistics file has incompatible columns: %s', ...
        statsFile);

    nMarkers = 0;
    return;
end

% Current axis range
yl = ylim(ax);
yRange = diff(yl);

if ~isfinite(yRange) || yRange <= 0
    yRange = 1;
end

nMarkers = 0;
maxPlottedY = yl(2);

preColor = [0.72 0.62 0.85];
postColor = [0.36 0.18 0.55];

% Draw one chance line
locationText = lower(strtrim(string(T.Location)));
timeText = lower(strtrim(string(T.Time)));

locationMasks = { ...
    locationText == "in" | locationText == "inside" | ...
        locationText == "inside tank", ...
    locationText == "out" | locationText == "outside" | ...
        locationText == "outside tank"};

xRanges = [0.55 3.45; ... % Inside-tank region
           4.55 7.45];    % Outside-tank region

phaseMasks = { ...
    timeText == "pre" | contains(timeText,"pre-switch"), ...
    timeText == "post" | contains(timeText,"post-switch")};

phaseColors = {preColor,postColor};
phaseStyles = {':','--'}; % Pre = dotted; Post = dashed

for locationIndex = 1:2
    for phaseIndex = 1:2
        selectedRows = locationMasks{locationIndex} & ...
            phaseMasks{phaseIndex};
        chanceValues = double(T.ChanceMean(selectedRows));
        chanceValues = chanceValues(isfinite(chanceValues));

        if isempty(chanceValues)
            warning('No chance mean found for location %d, phase %d.', ...
                locationIndex,phaseIndex);
            continue;
        end

        chanceY = mean(chanceValues);
        plot(ax,xRanges(locationIndex,:),[chanceY chanceY], ...
            'Color',phaseColors{phaseIndex}, ...
            'LineStyle',phaseStyles{phaseIndex}, ...
            'LineWidth',1.5, ...
            'HandleVisibility','off');

        maxPlottedY = max(maxPlottedY,chanceY);
    end
end

for row = 1:height(T)

    % Read condition information
    location = lower(strtrim(string(T.Location(row))));
    shape = lower(strtrim(string(T.RobotShape(row))));
    time = lower(strtrim(string(T.Time(row))));

    % Map presentation location
    if location == "in" || contains(location,"inside")

        locationOffset = 0;

    elseif location == "out" || contains(location,"outside")

        locationOffset = 3;

    else

        warning( ...
            'Could not map chance-comparison location: %s', ...
            location);

        continue;
    end

    % Map stimulus morphology
    if shape == "rod"

        conditionIndex = locationOffset + 1;

    elseif shape == "rectangle"

        conditionIndex = locationOffset + 2;

    elseif shape == "fish" || ...
           shape == "fish model" || ...
           shape == "fishmodel" || ...
           shape == "3d replica"

        conditionIndex = locationOffset + 3;

    else

        warning( ...
            'Could not map stimulus morphology: %s', ...
            shape);

        continue;
    end

    % Map time phase
    if time == "pre" || contains(time,"pre-switch")

        phaseIndex = 1;

    elseif time == "post" || contains(time,"post-switch")

        phaseIndex = 2;

    else

        warning( ...
            'Could not map time phase: %s', ...
            time);

        continue;
    end

    markerX = ...
        barHandles(phaseIndex).XEndPoints(conditionIndex);

    pFdr = double(T.p_FDR(row));

    if ~isfinite(pFdr) || pFdr >= 0.05
        continue;
    end

    markerY = ...
        meanValues(conditionIndex,phaseIndex) + ...
        semValues(conditionIndex,phaseIndex) + ...
        0.025*yRange;

    text( ...
        ax, ...
        markerX, ...
        markerY, ...
        '#', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',13, ...
        'FontWeight','bold', ...
        'Interpreter','none', ...
        'HandleVisibility','off');

    maxPlottedY = max(maxPlottedY,markerY);
    nMarkers = nMarkers + 1;

end

if maxPlottedY > yl(2)

    ylim( ...
        ax, ...
        [yl(1), maxPlottedY + 0.05*yRange]);

end

end
