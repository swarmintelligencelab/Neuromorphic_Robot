%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FishPool,RobotPool] = ...
    pool_conditions_for_interaction( ...
        TimeSeries_Fish,TimeSeries_Robot,idx)

% Determine total number of trials
nTotal = 0;
maxRows = 0;

for k = 1:numel(idx)

    c = idx(k);

    nTotal = nTotal + size(TimeSeries_Fish{c}.XX_f,2);

    maxRows = max( ...
        maxRows, ...
        max([ ...
            size(TimeSeries_Fish{c}.XX_f,1), ...
            size(TimeSeries_Robot{c}.XX_r,1) ...
        ]) ...
    );
end

% Allocate pooled matrices
FishPool.XX_f = nan(maxRows,nTotal);
FishPool.YY_f = nan(maxRows,nTotal);

RobotPool.XX_r = nan(maxRows,nTotal);
RobotPool.YY_r = nan(maxRows,nTotal);

RobotPool.SwitchIdx = nan(1,nTotal);
RobotPool.NumFrames = nan(1,nTotal);

% Pool trials
colStart = 1;

for k = 1:numel(idx)

    c = idx(k);

    Fish = TimeSeries_Fish{c};
    Robot = TimeSeries_Robot{c};

    nTrials = size(Fish.XX_f,2);

    if size(Robot.XX_r,2) ~= nTrials
        error( ...
            'Fish/robot trial mismatch in condition %d.', ...
            c);
    end

    cols = colStart:(colStart+nTrials-1);

    Tf = size(Fish.XX_f,1);
    Tr = size(Robot.XX_r,1);

    FishPool.XX_f(1:Tf,cols) = Fish.XX_f;
    FishPool.YY_f(1:Tf,cols) = Fish.YY_f;

    RobotPool.XX_r(1:Tr,cols) = Robot.XX_r;
    RobotPool.YY_r(1:Tr,cols) = Robot.YY_r;

    RobotPool.SwitchIdx(cols) = ...
        reshape(Robot.SwitchIdx,1,[]);

    RobotPool.NumFrames(cols) = ...
        reshape(Robot.NumFrames,1,[]);

    colStart = colStart + nTrials;

end

end
