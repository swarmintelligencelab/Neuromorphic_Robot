%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/06/2026
%Description : Chance interaction analysis using shuffled experimental fish-robot pairs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [chance_freq, chance_meanDur, chance_totalTime, ...
          freq_all, mean_all, total_all] = ...
    ComputeChanceFre_Dur_T_ShuffledGroup(XX_f, YY_f, XX_r, YY_r, dt, R_threshold, M)

[Tf,Nf] = size(XX_f);
[Tr,Nr] = size(XX_r);

assert(Nf == Nr,'Number of fish and robot trials must be equal.');
assert(Nf > 1,'At least two trials are required.');

N = Nf;

T = min([Tf,Tr,size(YY_f,1),size(YY_r,1)]);

XX_f = XX_f(1:T,:);
YY_f = YY_f(1:T,:);
XX_r = XX_r(1:T,:);
YY_r = YY_r(1:T,:);

% Define Pre and Post phases
split_idx = round(5*60/dt);
split_idx = min(split_idx,T);

phase_idx = {1:split_idx,(split_idx+1):T};

% Allocate null distributions
% Rows    = permutations
% Columns = Pre, Post
freq_all  = nan(M,2);
mean_all  = nan(M,2);
total_all = nan(M,2);

% Permutations
for m = 1:M

    robot_perm = randperm(N);

    while any(robot_perm == 1:N)
        robot_perm = randperm(N);
    end

    % Analyze Pre and Post separately
    for ph = 1:2

        idx = phase_idx{ph};

        if isempty(idx)
            continue
        end

        % One value for each shuffled fish-robot pair
        pair_freq      = nan(N,1);
        pair_meanDur   = nan(N,1);
        pair_totalTime = nan(N,1);


        % Analyze all fish
        for i = 1:N

            % Robot assigned to fish i in this permutation
            j = robot_perm(i);

            % Fish trajectory from trial i
            xf = XX_f(idx,i);
            yf = YY_f(idx,i);

            % Robot trajectory from a different trial j
            xr = XX_r(idx,j);
            yr = YY_r(idx,j);

            % Fish-robot distance
            dist = hypot(xf-xr,yf-yr);

            is_near = (dist <= R_threshold);

            is_near = is_near(:);

            % Detect interaction episodes
            visit_starts = find(diff([0;is_near]) == 1);

            visit_ends = find(diff([is_near;0]) == -1);

            if ~isempty(visit_starts) && ~isempty(visit_ends)

                if visit_ends(1) < visit_starts(1)
                    visit_ends(1) = [];
                end

                if length(visit_starts) > length(visit_ends)
                    visit_starts(end) = [];
                end

            end

            % Episode durations
            visit_durations = (visit_ends-visit_starts+1)*dt;

            n_ep = numel(visit_durations);

            Tsec = numel(idx)*dt;

            % Interaction frequency
            pair_freq(i) = n_ep/Tsec;

            % Mean interaction duration
            % If no interaction occurred, duration is undefined -> NaN
            if n_ep > 0

                pair_meanDur(i) = mean(visit_durations,'omitnan');

            else

                pair_meanDur(i) = NaN;

            end

            % Total interaction time
            % Zero is valid if no interaction occurred.
            pair_totalTime(i) = sum(is_near)*dt;

        end

        % Group mean for this permutation
        freq_all(m,ph)  = mean(pair_freq,'omitnan');

        mean_all(m,ph)  = mean(pair_meanDur,'omitnan');

        total_all(m,ph) = mean(pair_totalTime,'omitnan');

    end

end

chance_freq = mean(freq_all(:),'omitnan');
chance_meanDur = mean(mean_all(:),'omitnan');
chance_totalTime = mean(total_all(:),'omitnan');

end
