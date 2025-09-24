%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [freq_matrix_C2A, mean_dur_matrix_C2A, total_time_matrix_C2A, ...
          freq_matrix_A2C, mean_dur_matrix_A2C, total_time_matrix_A2C] = ...
    FrequencyDuration_InspectionEpisodes(TimeSeries_Fish_C2A, TimeSeries_Fish_A2C, ...
                                         TimeSeries_Robot_C2A, TimeSeries_Robot_A2C, ...
                                         dt, BL, proximity_factor)

R_threshold = proximity_factor * BL;
split_idx = round(5 * 60 / dt);  % 5 minutes

% Container for outputs
freq_matrix_C2A = [];
mean_dur_matrix_C2A = [];
total_time_matrix_C2A = [];

freq_matrix_A2C = [];
mean_dur_matrix_A2C = [];
total_time_matrix_A2C = [];

conditions = {'C2A', 'A2C'};
fish_data = {TimeSeries_Fish_C2A, TimeSeries_Fish_A2C};
robot_data = {TimeSeries_Robot_C2A, TimeSeries_Robot_A2C};

for condIdx = 1:2
    Fish = fish_data{condIdx};
    Robot = robot_data{condIdx};

    XX_f = Fish.XX_f;
    YY_f = Fish.YY_f;
    XX_r = Robot.XX_r;
    YY_r = Robot.YY_r;

    n_trials = size(XX_f, 2);

    freq_mat = NaN(2, n_trials);
    dur_mat = NaN(2, n_trials);
    time_mat = NaN(2, n_trials);

    for trial = 1:n_trials
        for phase = 1:2  % 1 = Before, 2 = After
            if phase == 1
                idx = 1:split_idx;
            else
                idx = (split_idx+1):length(XX_f(:,trial));
            end

            xf = XX_f(idx, trial); yf = YY_f(idx, trial);
            xr = XX_r(idx, trial); yr = YY_r(idx, trial);

            dist = sqrt((xf - xr).^2 + (yf - yr).^2);
            is_near = dist <= R_threshold;

            visit_starts = find(diff([0; is_near]) == 1);
            visit_ends   = find(diff([is_near; 0]) == -1);
            if ~isempty(visit_starts) && ~isempty(visit_ends)
                if visit_ends(1) < visit_starts(1)
                    visit_ends(1) = [];
                end
                if length(visit_starts) > length(visit_ends)
                    visit_starts(end) = [];
                end
            end

            visit_durations = (visit_ends - visit_starts) * dt;

            freq_mat(phase, trial) = numel(visit_starts) / (length(is_near) * dt);
            dur_mat(phase, trial) = mean(visit_durations);
            time_mat(phase, trial) = sum(is_near) * dt;
        end
    end

    % Assign to output variables
    if condIdx == 1
        freq_matrix_C2A = freq_mat;
        mean_dur_matrix_C2A = dur_mat;
        total_time_matrix_C2A = time_mat;
    else
        freq_matrix_A2C = freq_mat;
        mean_dur_matrix_A2C = dur_mat;
        total_time_matrix_A2C = time_mat;
    end
    
end