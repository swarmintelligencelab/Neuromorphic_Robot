%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [chance_freq, chance_meanDur, chance_totalTime, ...
          freq_all, mean_all, total_all] = ...
    compute_chance_interactions_shuffled_group( ...
        Fish, Robot, Parameters, R_threshold, M)

XX_f = Fish.XX_f;
YY_f = Fish.YY_f;
XX_r = Robot.XX_r;
YY_r = Robot.YY_r;

nTrials = size(XX_f,2);

if size(XX_r,2) ~= nTrials
    error('Fish and robot must have the same number of trials.');
end

freq_all  = nan(M,2);
mean_all  = nan(M,2);
total_all = nan(M,2);

for m = 1:M

    robotOrder = randperm(nTrials);

    while any(robotOrder == 1:nTrials)
        robotOrder = randperm(nTrials);
    end

    % One row per shuffled fish-robot pair
    freq_pairs  = nan(nTrials,2);
    mean_pairs  = nan(nTrials,2);
    total_pairs = nan(nTrials,2);

    for i = 1:nTrials

        j = robotOrder(i);

        % Fish trial i uses the switch timing from its original trial
        fishPhase = switch_phase_indices( ...
            Robot.SwitchIdx(i), ...
            Robot.NumFrames(i), ...
            Parameters);

        % Shuffled robot trial j uses its own switch timing
        robotPhase = switch_phase_indices( ...
            Robot.SwitchIdx(j), ...
            Robot.NumFrames(j), ...
            Parameters);

        for phase = 1:2

            idxF = fishPhase{phase};
            idxR = robotPhase{phase};

            if isempty(idxF) || isempty(idxR)
                continue;
            end

            % Align the two phase segments
            n = min(numel(idxF),numel(idxR));

            idxF = idxF(1:n);
            idxR = idxR(1:n);

            xf = XX_f(idxF,i);
            yf = YY_f(idxF,i);

            xr = XX_r(idxR,j);
            yr = YY_r(idxR,j);

            valid = ~isnan(xf) & ~isnan(yf) & ...
                    ~isnan(xr) & ~isnan(yr);

            xf = xf(valid);
            yf = yf(valid);
            xr = xr(valid);
            yr = yr(valid);

            if isempty(xf)
                continue;
            end

            dist = hypot(xf-xr, yf-yr);
            is_near = dist <= R_threshold;

            [visit_starts, visit_ends] = episode_edges(is_near);

            visit_durations = ...
                (visit_ends-visit_starts+1) * Parameters.dt;

            observationTime = numel(is_near) * Parameters.dt;

            % Frequency
            freq_pairs(i,phase) = ...
                numel(visit_starts) / observationTime;

            % Mean duration
            if isempty(visit_durations)
                mean_pairs(i,phase) = 0;
            else
                mean_pairs(i,phase) = ...
                    mean(visit_durations,'omitnan');
            end

            % Total interaction time
            total_pairs(i,phase) = ...
                sum(is_near) * Parameters.dt;

        end
    end

    % Mean over all shuffled pairs in this permutation
    freq_all(m,:) = mean(freq_pairs,1,'omitnan');
    mean_all(m,:) = mean(mean_pairs,1,'omitnan');
    total_all(m,:) = mean(total_pairs,1,'omitnan');

end

% Mean of null distributions = chance level
chance_freq = mean(freq_all,1,'omitnan');
chance_meanDur = mean(mean_all,1,'omitnan');
chance_totalTime = mean(total_all,1,'omitnan');

end
