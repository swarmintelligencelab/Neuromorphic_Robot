%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/06/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Mean_matrix] = Compute_time_Average(M,num_bins)

N = size(M, 1); 
num_individuals = size(M, 2); % Number of individuals

% Compute bin size
bin_size = N / num_bins; % Samples per bin

% Initialize storage
mean_vals = zeros(num_bins, 1);
var_vals = zeros(num_bins, 1);
sem_vals = zeros(num_bins, 1); % Standard Error of the Mean

% Compute statistics per bin
for i = 1:num_bins
    idx_start = (i-1) * bin_size + 1;
    idx_end = i * bin_size;
    
    % Extract data for the current bin
    bin_data = M(idx_start:idx_end, :); % Bin range for all individuals
    
    % Compute mean and variance across individuals for this bin
    Mean_matrix(i,:) = nanmean(bin_data, 1);
end

end
