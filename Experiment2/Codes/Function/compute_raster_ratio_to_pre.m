%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ratioByCondition = compute_raster_ratio_to_pre(rawByCondition, plotTimeMin)

preMask = plotTimeMin < 5;
ratioByCondition = cell(size(rawByCondition));

for c = 1:numel(rawByCondition)
    M = rawByCondition{c};
    preMean = mean(M(:, preMask), 2, 'omitnan');
    preMean(abs(preMean) < eps) = NaN;
    ratioByCondition{c} = M ./ preMean;
end

end
