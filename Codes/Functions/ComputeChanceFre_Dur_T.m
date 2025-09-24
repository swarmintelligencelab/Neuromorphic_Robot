%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [chance_freq, chance_meanDur, chance_totalTime] = ...
    ComputeChanceFre_Dur_T(XX_f, YY_f, XX_r, YY_r, dt, R_threshold, M)

[Tf, Nf] = size(XX_f);
[Tr, Nr] = size(XX_r);
T = min([Tf, Tr, size(YY_f,1), size(YY_r,1)]);
XX_f = XX_f(1:T,:); YY_f = YY_f(1:T,:);
XX_r = XX_r(1:T,:); YY_r = YY_r(1:T,:);

split_idx = round(5 * 60 / dt);
split_idx = min(split_idx, T);
phase_idx = {1:split_idx, (split_idx+1):T};

freq_all  = zeros(M,1);
mean_all  = zeros(M,1);
total_all = zeros(M,1);

for m = 1:M
    % Random select one fish and one robot
    i = randi(Nf); j = randi(Nr);
    xf_full = XX_f(:,i); yf_full = YY_f(:,i);
    xr_full = XX_r(:,j); yr_full = YY_r(:,j);

    freq_phase = nan(2,1);
    mean_phase = nan(2,1);
    time_phase = nan(2,1);

    for ph = 1:2
        idx = phase_idx{ph};
        if isempty(idx), continue; end

        xf = xf_full(idx); yf = yf_full(idx);
        xr = xr_full(idx); yr = yr_full(idx);

        dist = hypot(xf - xr, yf - yr);
        is_near = (dist <= R_threshold);

        visit_starts = find(diff([0; is_near]) == 1);
        visit_ends   = find(diff([is_near; 0]) == -1);
        if ~isempty(visit_starts) && ~isempty(visit_ends)
            if visit_ends(1) < visit_starts(1)
                visit_ends(1) = [];
            end
            if length(visit_starts) > length(visit_ends)
                visit_starts(end) = [];
            end
        end

        % Duration
        visit_durations = (visit_ends - visit_starts) * dt;

        n_ep  = numel(visit_durations);
        Tsec  = numel(idx) * dt;

        freq_phase(ph) = n_ep / Tsec;
        mean_phase(ph) = mean(visit_durations,'omitnan');
        time_phase(ph) = sum(is_near) * dt;
    end

    freq_all(m)  = mean(freq_phase, 'omitnan');
    mean_all(m)  = mean(mean_phase, 'omitnan');
    total_all(m) = mean(time_phase, 'omitnan');
end

chance_freq      = mean(freq_all, 'omitnan');
chance_meanDur   = mean(mean_all, 'omitnan');
chance_totalTime = mean(total_all, 'omitnan');

end