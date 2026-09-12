%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function idx = switch_phase_indices(switchIdx, numFrames, Parameters)

idx = {[], []};

if ~isfinite(switchIdx) || ~isfinite(numFrames)
    return;
end

switchIdx = round(switchIdx);
numFrames = round(numFrames);
preFrames = Parameters.phaseWindowFrames;
postFrames = Parameters.phaseWindowFrames;

preStart = switchIdx - preFrames;
preEnd = switchIdx - 1;
postStart = switchIdx;
postEnd = switchIdx + postFrames - 1;

hasFullPre = preStart >= 1;
hasFullPost = postEnd <= numFrames;

if Parameters.requireFullPhaseWindow
    if hasFullPre
        idx{1} = preStart:preEnd;
    end
    if hasFullPost
        idx{2} = postStart:postEnd;
    end
else
    idx{1} = max(1, preStart):min(numFrames, preEnd);
    idx{2} = max(1, postStart):min(numFrames, postEnd);
end

end
