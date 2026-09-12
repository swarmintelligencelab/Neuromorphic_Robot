%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function order = grouped_condition_plot_order(conditionDefs)

objectOrder = {'Rod', 'Rectangle', 'FishModel'};
locationOrder = {'In', 'Out'};
order = [];

for iObj = 1:numel(objectOrder)
    for iLoc = 1:numel(locationOrder)
        idx = find(strcmp({conditionDefs.object}, objectOrder{iObj}) & ...
            strcmp({conditionDefs.location}, locationOrder{iLoc}), 1);
        if ~isempty(idx)
            order(end + 1) = idx;
        end
    end
end

end
