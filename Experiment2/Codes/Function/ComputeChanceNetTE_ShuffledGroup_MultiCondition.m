%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [chance_TE, TE_null, TAU_null] = ...
    ComputeChanceNetTE_ShuffledGroup_MultiCondition( ...
        B_f, B_r, nbin, delay_values, M)

[nLocations,nPhases] = size(B_f);

if nPhases ~= 2
    error('B_f and B_r must have two phase columns: Pre and Post.');
end

if ~isequal(size(B_f),size(B_r))
    error('B_f and B_r must have the same cell-array size.');
end

TE_null = cell(nLocations,nPhases);

for g = 1:nLocations
    for phase = 1:nPhases
        TE_null{g,phase} = nan(M,1);
    end
end

TAU_null = nan(M,1);


for m = 1:M

    B_r_shuffle = cell(nLocations,nPhases);

    for g = 1:nLocations

        nTrials = size(B_f{g,1},2);

        if nTrials < 2
            error( ...
                'Location %d needs at least two trials for shuffling.', ...
                g);
        end

        for phase = 1:nPhases

            if size(B_f{g,phase},2) ~= nTrials || ...
               size(B_r{g,phase},2) ~= nTrials
                error( ...
                    'Pre/Post fish/robot trial counts differ in location %d.', ...
                    g);
            end
        end

        robotOrder = randperm(nTrials);

        while any(robotOrder == 1:nTrials)
            robotOrder = randperm(nTrials);
        end

        for phase = 1:nPhases
            B_r_shuffle{g,phase} = ...
                B_r{g,phase}(:,robotOrder);
        end
    end

    % Evaluate every candidate TAU
    TE_tau = nan(1,numel(delay_values));

    TE_group_tau = cell( ...
        nLocations, ...
        nPhases, ...
        numel(delay_values));

    for k = 1:numel(delay_values)

        tauNow = delay_values(k);
        TE_all = [];

        for g = 1:nLocations
            for phase = 1:nPhases

                tmp = NetTransferEntropy_main( ...
                    B_f{g,phase}, ...
                    B_r_shuffle{g,phase}, ...
                    nbin, ...
                    tauNow);

                TE_group_tau{g,phase,k} = tmp;
                TE_all = [TE_all; tmp(:)];
            end
        end

        TE_tau(k) = mean(TE_all,'omitnan');
    end

    % Common optimal TAU for this permutation
    [~,tauIndex] = max(TE_tau);

    TAU_null(m) = delay_values(tauIndex);

    % Store one null group mean for every Location x Phase
    for g = 1:nLocations
        for phase = 1:nPhases

            tmp = TE_group_tau{g,phase,tauIndex};

            TE_null{g,phase}(m) = ...
                mean(tmp,'omitnan');
        end
    end
end

% Mean chance level for every Location x Phase
chance_TE = cell(nLocations,nPhases);

for g = 1:nLocations
    for phase = 1:nPhases
        chance_TE{g,phase} = ...
            mean(TE_null{g,phase},'omitnan');
    end
end

end
