%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function metricCells = build_te_metric_cells(TE,conditionDefs)

requiredFields = { ...
    'PoolNames', ...
    'NetTEByLocationPhase', ...
    'ConditionByLocation'};

assert(all(isfield(TE,requiredFields)), ...
    'Results.TE is missing fields required for the TE barplot.');

poolNames = string(TE.PoolNames);
nConditions = numel(conditionDefs);
metricCells = cell(1,nConditions);

for c = 1:nConditions

    conditionName = string(conditionDefs(c).display);

    if startsWith(conditionName,"In-",'IgnoreCase',true)
        locationName = "In";
    elseif startsWith(conditionName,"Out-",'IgnoreCase',true)
        locationName = "Out";
    else
        error('Cannot determine TE location for condition %s.',conditionName);
    end

    poolIndex = find(strcmpi(poolNames,locationName),1);
    assert(~isempty(poolIndex), ...
        'TE pool %s was not found in Results.TE.PoolNames.',locationName);

    conditionMetadata = string(TE.ConditionByLocation{poolIndex}(:));
    preValues = TE.NetTEByLocationPhase{poolIndex,1}(:);
    postValues = TE.NetTEByLocationPhase{poolIndex,2}(:);

    assert(numel(preValues) == numel(conditionMetadata) && ...
           numel(postValues) == numel(conditionMetadata), ...
        'TE values and condition metadata have different lengths for %s.', ...
        locationName);

    conditionMask = strcmpi(conditionMetadata,conditionName);
    assert(any(conditionMask), ...
        'No trial-level TE values were found for condition %s.',conditionName);

    metricCells{c} = [ ...
        preValues(conditionMask).'; ...
        postValues(conditionMask).'];
end

end
