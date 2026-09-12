%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Mean_matrix = compute_switch_phase_average(M, phaseInfo, Parameters)

nTrials = size(M, 2);
Mean_matrix = nan(2, nTrials);

for trial = 1:nTrials
    idx = switch_phase_indices(phaseInfo.SwitchIdx(trial), ...
        phaseInfo.NumFrames(trial), Parameters);

    for phase = 1:2
        if isempty(idx{phase}), continue; end
        Mean_matrix(phase, trial) = mean(M(idx{phase}, trial), 'omitnan');
    end
end

end
