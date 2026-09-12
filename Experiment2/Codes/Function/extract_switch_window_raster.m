%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function M = extract_switch_window_raster(signalMatrix, windowFrames)

nTrials = size(signalMatrix, 2);
M = nan(nTrials, windowFrames);
writeFrames = min(windowFrames, size(signalMatrix, 1));
M(:, 1:writeFrames) = signalMatrix(1:writeFrames, :)';

end
