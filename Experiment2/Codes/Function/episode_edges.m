%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [visit_starts, visit_ends] = episode_edges(is_near)

is_near = is_near(:);
visit_starts = find(diff([0; is_near]) == 1);
visit_ends = find(diff([is_near; 0]) == -1);

if ~isempty(visit_starts) && ~isempty(visit_ends)
    if visit_ends(1) < visit_starts(1)
        visit_ends(1) = [];
    end
    if numel(visit_starts) > numel(visit_ends)
        visit_starts(end) = [];
    end
end

end
