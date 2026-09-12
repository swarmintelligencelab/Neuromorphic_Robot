%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function assert_required_data(S, conditionDefs)

for c = 1:numel(conditionDefs)
    datName = ['Dat_' conditionDefs(c).name];
    if ~isfield(S, datName)
        error('Data file is missing required variable: %s', datName);
    end
end

end
