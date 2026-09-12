%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Bin_matrix = compute_switch_timebin_average(M, phaseInfo, Parameters, numBins)

nTrials = size(M, 2);
Bin_matrix = nan(numBins, nTrials);

if mod(numBins, 2) ~= 0
    error('Switch-aligned time bins require an even number of bins.');
end

binsPerPhase = numBins / 2;
framesPerBin = floor(Parameters.phaseWindowFrames / binsPerPhase);

for trial = 1:nTrials
    idx = switch_phase_indices(phaseInfo.SwitchIdx(trial), ...
        phaseInfo.NumFrames(trial), Parameters);

    if isempty(idx{1}) || isempty(idx{2})
        continue;
    end

    preIdx = idx{1};
    postIdx = idx{2};

    for b = 1:binsPerPhase
        localIdx = (b - 1) * framesPerBin + (1:framesPerBin);
        Bin_matrix(b, trial) = mean(M(preIdx(localIdx), trial), 'omitnan');
        Bin_matrix(b + binsPerPhase, trial) = mean(M(postIdx(localIdx), trial), 'omitnan');
    end
end

end
