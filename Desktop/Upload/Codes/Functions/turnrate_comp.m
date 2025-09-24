%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [WG] = turnrate_comp(costheta, sintheta, dt)

posonCircle = [costheta, sintheta];
WG = zeros(length(costheta)-1,1);
for k = 1 : length(costheta)-1
    v1 = posonCircle(k,:);
    v2 = posonCircle(k+1,:);
    deltaPhi = atan2(det([v1',v2']), dot(v1,v2));
    deltaPhi = wrapToPi(deltaPhi);
    WG(k) = (deltaPhi / dt);
end

end
