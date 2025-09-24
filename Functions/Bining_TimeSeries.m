%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%Description : Bin time‐series data into 2–5 discrete states based on standard-deviation thresholds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Q_f, Q_r] = Bining_TimeSeries(AA_f, AA_r, e_max, n_bins)

% Validate input
if n_bins < 2 || n_bins > 5
    error('n_bins must be an integer between 2 and 5');
end

% Compute base threshold using standard deviation
epsilon = e_max * std(AA_f(:), 'omitnan');

% Initialize output
Q_f = zeros(size(AA_f));
Q_r = zeros(size(AA_r));

% Apply binning
switch n_bins
    case 2
        % Binary: <0 and ≥0
        Q_f(AA_f < 0) = 1;
        Q_f(AA_f >= 0) = 2;

        Q_r(AA_r < 0) = 1;
        Q_r(AA_r >= 0) = 2;

    case 3
        % Classic: decel, dead zone, accel
        Q_f(AA_f < -epsilon) = 1;
        Q_f(abs(AA_f) <= epsilon) = 2;
        Q_f(AA_f > epsilon) = 3;

        Q_r(AA_r < -epsilon) = 1;
        Q_r(abs(AA_r) <= epsilon) = 2;
        Q_r(AA_r > epsilon) = 3;

    case 4
        % Add weak decel and weak accel
        Q_f(AA_f < -epsilon) = 1;
        Q_f(AA_f >= -epsilon & AA_f < 0) = 2;
        Q_f(AA_f >= 0 & AA_f <= epsilon) = 3;
        Q_f(AA_f > epsilon) = 4;

        Q_r(AA_r < -epsilon) = 1;
        Q_r(AA_r >= -epsilon & AA_r < 0) = 2;
        Q_r(AA_r >= 0 & AA_r <= epsilon) = 3;
        Q_r(AA_r > epsilon) = 4;

    case 5
        % Strong/weak decel, dead zone, weak/strong accel
        Q_f(AA_f < -2*epsilon) = 1;
        Q_f(AA_f >= -2*epsilon & AA_f < -epsilon) = 2;
        Q_f(abs(AA_f) <= epsilon) = 3;
        Q_f(AA_f > epsilon & AA_f <= 2*epsilon) = 4;
        Q_f(AA_f > 2*epsilon) = 5;

        Q_r(AA_r < -2*epsilon) = 1;
        Q_r(AA_r >= -2*epsilon & AA_r < -epsilon) = 2;
        Q_r(abs(AA_r) <= epsilon) = 3;
        Q_r(AA_r > epsilon & AA_r <= 2*epsilon) = 4;
        Q_r(AA_r > 2*epsilon) = 5;
end

end