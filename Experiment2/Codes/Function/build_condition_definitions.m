%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function conditionDefs = build_condition_definitions()

conditionDefs = struct( ...
    'name', {}, 'location', {}, 'object', {}, 'display', {}, 'color', {});

conditionDefs(1).name = 'in_robot';
conditionDefs(1).location = 'In';
conditionDefs(1).object = 'FishModel';
conditionDefs(1).display = 'In-FishModel';
conditionDefs(1).color = [1.00 0.40 0.40];

conditionDefs(2).name = 'in_rec';
conditionDefs(2).location = 'In';
conditionDefs(2).object = 'Rectangle';
conditionDefs(2).display = 'In-Rectangle';
conditionDefs(2).color = [0.20 0.60 1.00];

conditionDefs(3).name = 'in_rod';
conditionDefs(3).location = 'In';
conditionDefs(3).object = 'Rod';
conditionDefs(3).display = 'In-Rod';
conditionDefs(3).color = [0.55 0.55 0.55];

conditionDefs(4).name = 'out_robot';
conditionDefs(4).location = 'Out';
conditionDefs(4).object = 'FishModel';
conditionDefs(4).display = 'Out-FishModel';
conditionDefs(4).color = [0.75 0.10 0.10];

conditionDefs(5).name = 'out_rec';
conditionDefs(5).location = 'Out';
conditionDefs(5).object = 'Rectangle';
conditionDefs(5).display = 'Out-Rectangle';
conditionDefs(5).color = [0.00 0.25 0.75];

conditionDefs(6).name = 'out_rod';
conditionDefs(6).location = 'Out';
conditionDefs(6).object = 'Rod';
conditionDefs(6).display = 'Out-Rod';
conditionDefs(6).color = [0.05 0.05 0.05];

end
