%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Parameters = load_constants_6conditions(S, Metadata, conditionDefs)

x_max = 76;    % cm
y_tank = 50.4; % cm
y_max = 41;    % cm, bottom to water height

if isfield(Metadata, 'fps')
    fps = Metadata.fps;
else
    fps = 30;
end

if isfield(Metadata, 'dt')
    dt = Metadata.dt;
else
    dt = 1 / fps;
end

if isfield(Metadata, 'analysisFrames')
    T_total = Metadata.analysisFrames;
else
    firstDat = S.(['Dat_' conditionDefs(1).name]);
    T_total = numel(firstDat(1).Xf);
end

% Calibrate pixels to tank-centered cm from pooled fish trajectories.
minX = nan(1, numel(conditionDefs));
maxX = nan(1, numel(conditionDefs));
minY = nan(1, numel(conditionDefs));
maxY = nan(1, numel(conditionDefs));

for c = 1:numel(conditionDefs)
    Dat = S.(['Dat_' conditionDefs(c).name]);
    allX = vertcat_column(Dat, 'Xf');
    allY = vertcat_column(Dat, 'Yf');
    minX(c) = min(allX, [], 'omitnan');
    maxX(c) = max(allX, [], 'omitnan');
    minY(c) = min(allY, [], 'omitnan');
    maxY(c) = max(allY, [], 'omitnan');
end

Xmin = median(minX, 'omitnan');
Xmax = median(maxX, 'omitnan');
Ymin = median(minY, 'omitnan');
Ymax = median(maxY, 'omitnan');

scale_x = x_max / 2.1;
scale_y = y_max / 2.1;

x_rect = [-x_max/2, x_max/2, x_max/2, -x_max/2, -x_max/2];
y_rect = [-y_max/2, -y_max/2, y_max/2, y_max/2, -y_max/2];

Parameters.x_max = x_max;
Parameters.y_tank = y_tank;
Parameters.y_max = y_max;
Parameters.fps = fps;
Parameters.dt = dt;
Parameters.Xmin = Xmin;
Parameters.Xmax = Xmax;
Parameters.Ymin = Ymin;
Parameters.Ymax = Ymax;
Parameters.scale_x = scale_x;
Parameters.scale_y = scale_y;
Parameters.x_rect = x_rect;
Parameters.y_rect = y_rect;
Parameters.T_total = T_total;
Parameters.phaseWindowSec = 300;
Parameters.phaseWindowFrames = round(Parameters.phaseWindowSec / dt);
Parameters.requireFullPhaseWindow = true;

end
