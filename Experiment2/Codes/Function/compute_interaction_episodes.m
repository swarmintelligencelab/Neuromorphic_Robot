%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [freq_mat, dur_mat, time_mat, event_count_mat, observation_time_mat] = ...
    compute_interaction_episodes(Fish, Robot, Parameters, R_threshold)

XX_f = Fish.XX_f;
YY_f = Fish.YY_f;
XX_r = Robot.XX_r;
YY_r = Robot.YY_r;

nTrials = size(XX_f, 2);

freq_mat = nan(2, nTrials);
dur_mat = nan(2, nTrials);
time_mat = nan(2, nTrials);

event_count_mat = nan(2, nTrials);
observation_time_mat = nan(2, nTrials);

for trial = 1:nTrials

    phaseIdx = switch_phase_indices( ...
        Robot.SwitchIdx(trial), ...
        Robot.NumFrames(trial), ...
        Parameters);

    for phase = 1:2

        idx = phaseIdx{phase};

        if isempty(idx)
            continue;
        end

        xf = XX_f(idx, trial);
        yf = YY_f(idx, trial);
        xr = XX_r(idx, trial);
        yr = YY_r(idx, trial);

        dist = hypot(xf - xr, yf - yr);

        is_near = dist <= R_threshold;
        is_near(isnan(dist)) = false;

        % Detect interaction episodes
        [visit_starts, visit_ends] = episode_edges(is_near);

        visit_durations = ...
            (visit_ends - visit_starts + 1) * Parameters.dt;

        % Numerator: number of interaction episodes
        event_count_mat(phase, trial) = numel(visit_starts);

        % Denominator: total observation time in this phase
        observation_time_mat(phase, trial) = ...
            numel(is_near) * Parameters.dt;

        % Frequency
        freq_mat(phase, trial) = ...
            event_count_mat(phase, trial) / ...
            observation_time_mat(phase, trial);

        % Mean duration of interaction episodes
        if isempty(visit_durations)

            dur_mat(phase, trial) = 0;

        else

            dur_mat(phase, trial) = ...
                mean(visit_durations, 'omitnan');

        end

        % Total interaction time
        time_mat(phase, trial) = ...
            sum(is_near) * Parameters.dt;

    end
end

end
