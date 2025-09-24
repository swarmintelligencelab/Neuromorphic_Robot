%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function chance_TE = ComputeChanceTE(AA_f_Control, AA_f_C2A, AA_r_C2A, AA_f_A2C, AA_r_A2C, nbin, TAU, M)

    % Merge all individuals into a single dataset
    All_Data = [AA_f_Control, AA_f_C2A, AA_r_C2A, AA_f_A2C, AA_r_A2C];  
    [N, num_individuals] = size(All_Data); 
  
    % Initialize storage for chance TE values
    TE_chance_values = zeros(1, M);

    % Compute TE for M random interactions
    for m = 1:M
        % Randomly select two different individuals
        rand_inds = randperm(num_individuals, 2);
        B_X = All_Data(:, rand_inds(1)); % First random individual
        B_Y = All_Data(:, rand_inds(2)); % Second random individual
      
        % Compute TE
        TE_chance_values_a= TransEntropy(B_X, B_Y, nbin, TAU);
        TE_chance_values(m) = TE_chance_values_a;

    end

    % Compute mean and standard deviation of chance TE
    mean_chance_TE = mean(TE_chance_values);
    std_chance_TE = std(TE_chance_values);

    % Display results
    disp(['Mean Chance TE: ', num2str(mean_chance_TE)]);
    disp(['Std of Chance TE: ', num2str(std_chance_TE)]);

    % Return chance TE values
    chance_TE = TE_chance_values;
    
end