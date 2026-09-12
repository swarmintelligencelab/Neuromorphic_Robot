%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function climValues = compute_raster_robust_clim(dataByCondition, lowPct, highPct)

allValues = [];
for c = 1:numel(dataByCondition)
    M = dataByCondition{c};
    allValues = [allValues; M(:)];
end

allValues = allValues(isfinite(allValues));
if isempty(allValues)
    climValues = [0 1];
    return;
end

climValues = prctile(allValues, [lowPct, highPct]);
if diff(climValues) < 0.25
    centerValue = mean(climValues);
    climValues = centerValue + [-0.125 0.125];
end

if climValues(1) == climValues(2)
    climValues = climValues + [-0.5 0.5];
end

end
