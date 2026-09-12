%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/06/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [chance_TE, TE_null_all] = ComputeChanceNetTE_ShuffledGroup(B_f, B_r, nbin, TAU, M)

[Nt,Nf] = size(B_f);

[Ntr,Nr] = size(B_r);

assert(Nt == Ntr,'Fish and robot time series must have the same number of time samples.');

assert(Nf == Nr,'Number of fish and robot trials must be equal.');

assert(Nf > 1,'At least two trials are required.');

N = Nf;

TE_null_all = nan(M,1);

for m = 1:M

    robot_perm = randperm(N);

    while any(robot_perm == 1:N)

        robot_perm = randperm(N);

    end

    B_r_shuffled = B_r(:,robot_perm);

    NetTE_null = NetTransferEntropy_main(B_f,B_r_shuffled,nbin,TAU);

    TE_null_all(m) = mean(NetTE_null,'omitnan');

end

chance_TE = mean(TE_null_all,'omitnan');

end
