%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = TransEntropy(X, Y, bin, TAU)
    %% This function computes Transfer Entropy between two time series.
    % X, Y: Binned time series
    % bin: Number of bins
    % TAU: Time delay for TE computation

    SimulationSteps = length(X);

    % Initialize probability distribution
    CONT = zeros(bin, bin, bin);

    % Compute joint frequency counts
    for a = 1:bin
        for b = 1:bin
            for c = 1:bin
                for t = (TAU + 1):(SimulationSteps - 1)
                    if (t - TAU > 0) && (Y(t) == b) && (Y(t + 1) == a) && (X(t - TAU) == c)
                        CONT(a, b, c) = CONT(a, b, c) + 1;
                    end
                end
            end
        end
    end

    % Normalize probabilities
    CONT = CONT / sum(CONT(:));  % Ensures proper probability distribution

    % Compute Transfer Entropy
    TE = 0;
    for a = 1:bin
        for b = 1:bin
            for c = 1:bin
                num = CONT(a, b, c) * sum(sum(CONT(:, b, :)));
                den = sum(CONT(:, b, c)) * sum(CONT(a, b, :));

                if num > 0 && den > 0
                    TE = TE + CONT(a, b, c) * log2(num / den);
                end
            end
        end
    end

    y = TE;
end