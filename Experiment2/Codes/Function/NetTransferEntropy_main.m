%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function NetTE = NetTransferEntropy_main(B_f,B_r,nbin,TAU)

nFish = size(B_f,2);
nRobot = size(B_r,2);

if nFish ~= nRobot
    error( ...
        'B_f and B_r have different trial numbers: fish=%d, robot=%d.', ...
        nFish,nRobot);
end

nTrials = nFish;

NetTE = nan(1,nTrials);

for i = 1:nTrials

    X = B_r(:,i);
    Y = B_f(:,i);

    TE_rf = TransEntropy(X,Y,nbin,TAU);
    % TE_fr = TransEntropy(Y,X,nbin,TAU);

    NetTE(i) = TE_rf;

end

end
