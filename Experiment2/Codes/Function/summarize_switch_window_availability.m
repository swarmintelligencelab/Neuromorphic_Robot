%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function T = summarize_switch_window_availability(TimeSeries_Robot, TrialIDs, conditionDefs, Parameters)

rows = {};

for c = 1:numel(TimeSeries_Robot)
    Robot = TimeSeries_Robot{c};
    trialIDs = TrialIDs{c};

    for trial = 1:numel(trialIDs)
        switchIdx = Robot.SwitchIdx(trial);
        numFrames = Robot.NumFrames(trial);
        switchTime = Robot.SwitchTimeSec(trial);

        preSec = NaN;
        postSec = NaN;
        hasFullPre = false;
        hasFullPost = false;

        if isfinite(switchIdx) && isfinite(numFrames)
            preSec = (switchIdx - 1) * Parameters.dt;
            postSec = (numFrames - switchIdx + 1) * Parameters.dt;
            hasFullPre = preSec >= Parameters.phaseWindowSec;
            hasFullPost = postSec >= Parameters.phaseWindowSec;
        end

        rows(end+1, :) = { ...
            string(conditionDefs(c).name), ...
            string(conditionDefs(c).display), ...
            trialIDs(trial), ...
            switchTime, ...
            numFrames * Parameters.dt, ...
            preSec, ...
            postSec, ...
            hasFullPre, ...
            hasFullPost, ...
            hasFullPre && hasFullPost};
    end
end

T = cell2table(rows, 'VariableNames', { ...
    'ConditionName', 'ConditionLabel', 'TrialID', 'SwitchTimeSec', ...
    'UsedDurationSec', 'AvailablePreSec', 'AvailablePostSec', ...
    'HasFullPre300s', 'HasFullPost300s', 'HasFullPrePost300s'});

fprintf('\nSwitch-window availability, requiring %.0f s before and after switch:\n', ...
    Parameters.phaseWindowSec);

for c = 1:numel(conditionDefs)
    isCond = T.ConditionName == string(conditionDefs(c).name);
    nFull = nnz(T.HasFullPrePost300s(isCond));
    nTotal = nnz(isCond);
    fprintf('  %s: %d/%d trials have full pre/post windows\n', ...
        conditionDefs(c).display, nFull, nTotal);
end

end
