%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FishS, RobotS] = sample_pairs_from_control(FXX, FYY, RXX, RYY, n_trials, min_shift)

    T = size(FXX,1);
    fi = randi(size(FXX,2), 1, n_trials);   % Randomly select in Control
    ri = randi(size(RXX,2), 1, n_trials);   % Randomly select robots
    XX_f = zeros(T, n_trials); YY_f = zeros(T, n_trials);
    XX_r = zeros(T, n_trials); YY_r = zeros(T, n_trials);
    for k = 1:n_trials
        xf = FXX(:, fi(k));  yf = FYY(:, fi(k));
        xr = RXX(:, ri(k));  yr = RYY(:, ri(k));

        if min_shift > 0
            kshift = randi([min_shift, max(min_shift, T-1)]);
            xf = circshift(xf, kshift);
            yf = circshift(yf, kshift);
        end
        XX_f(:,k) = xf; YY_f(:,k) = yf;
        XX_r(:,k) = xr; YY_r(:,k) = yr;
    end
    FishS  = struct('XX_f', XX_f, 'YY_f', YY_f);
    RobotS = struct('XX_r', XX_r, 'YY_r', YY_r);

end