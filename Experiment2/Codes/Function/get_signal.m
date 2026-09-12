%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function signal = get_signal(DatTrial, fieldName, n)

signal = nan(n, 1);

if isfield(DatTrial, fieldName)
    raw = DatTrial.(fieldName)(:);
    write_len = min(n, numel(raw));
    signal(1:write_len) = raw(1:write_len);
end

end
