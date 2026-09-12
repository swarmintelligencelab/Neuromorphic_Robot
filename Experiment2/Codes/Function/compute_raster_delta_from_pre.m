%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function deltaByCondition = compute_raster_delta_from_pre(rawByCondition, plotTimeMin)

preMask = plotTimeMin < 5;
deltaByCondition = cell(size(rawByCondition));

for c = 1:numel(rawByCondition)
    M = rawByCondition{c};
    preMedian = median(M(:, preMask), 2, 'omitnan');
    deltaByCondition{c} = M - preMedian;
end

end
