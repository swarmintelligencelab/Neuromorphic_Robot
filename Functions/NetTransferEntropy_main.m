%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function NetTE = NetTransferEntropy_main(B_f,B_r,nbin,TAU)

% Compute Transfer Entropy for each individual
for i = 1:10
    X = B_r(:,i);
    Y = B_f(:,i);
    TE_rf = TransEntropy(X, Y, nbin, TAU);
    TE_fr = TransEntropy(Y, X, nbin, TAU);
    NetTE(i) = TE_rf - 0*TE_fr;
end

end