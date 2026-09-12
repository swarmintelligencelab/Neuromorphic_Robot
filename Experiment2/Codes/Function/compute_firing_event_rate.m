%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rate_mat,count_mat,exposure_mat] = ...
    compute_firing_event_rate(Robot,Parameters)

dt = Parameters.dt;
nTrials = size(Robot.NN_firing,2);

rate_mat     = nan(2,nTrials);
count_mat    = nan(2,nTrials);
exposure_mat = nan(2,nTrials);

binary_data = Robot.NN_firing == 1;

for tr = 1:nTrials

    nFrames = Robot.NumFrames(tr);

    x = binary_data(1:nFrames,tr);

    % Find firing-event onsets in the FULL trial
    event_onset = [x(1); diff(x) == 1];

    phaseIdx = switch_phase_indices( ...
        Robot.SwitchIdx(tr), ...
        nFrames, ...
        Parameters);

    for phase = 1:2

        idx = phaseIdx{phase};

        if isempty(idx)
            continue
        end

        count_mat(phase,tr) = sum(event_onset(idx));
        exposure_mat(phase,tr) = numel(idx)*dt;

        rate_mat(phase,tr) = ...
            count_mat(phase,tr) / exposure_mat(phase,tr);

    end
end

end
