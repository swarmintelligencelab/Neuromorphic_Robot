%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_grouped_2_bar(M_YY_f_A2C, M_YY_f_C2A, Font_Size)

    % Compute means and standard errors
    mean_A2C = mean(M_YY_f_A2C, 2);
    sem_A2C = std(M_YY_f_A2C, 0, 2) / sqrt(size(M_YY_f_A2C, 2));

    mean_C2A = mean(M_YY_f_C2A, 2);
    sem_C2A = std(M_YY_f_C2A, 0, 2) / sqrt(size(M_YY_f_C2A, 2));

    % Data arrangement for bar plot
    mean_values = [mean_A2C, mean_C2A]; 
    sem_values = [sem_A2C, sem_C2A]; 
    
    % Define colors for bars
    bar_colors = [1 0.4 0.4;     % Light Red for A2C
                  0.2 0.6 1];    % Blue for C2A
    
    % Create grouped bar plot
    figure;
    hold on;
    bar_handle = bar(mean_values, 'grouped'); 

    for k = 1:length(bar_handle)
        bar_handle(k).FaceColor = bar_colors(k, :);
        bar_handle(k).EdgeColor = 'none';
    end

    % Error bars
    num_groups = size(mean_values, 1); 
    num_bars = size(mean_values, 2); 

    x_positions = nan(num_groups, num_bars);
    for i = 1:num_bars
        x_positions(:, i) = bar_handle(i).XEndPoints;
    end

    errorbar(x_positions, mean_values, sem_values, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

    % Plot individual data points
    jitter_amount = 0.04;
    for g = 1:num_groups

        % A2C group
        x = x_positions(g,1) + (rand(1, size(M_YY_f_A2C, 2)) - 0.5) * jitter_amount;
        y = M_YY_f_A2C(g, :);

        % C2A group
        x = x_positions(g,2) + (rand(1, size(M_YY_f_C2A, 2)) - 0.5) * jitter_amount;
        y = M_YY_f_C2A(g, :);
    end

    % Formatting
    xticks(1:2);
    xticklabels({'Pre-Switch', ' Post-Switch'});
    set(gca, 'FontSize', Font_Size, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Color', 'w');
    grid off;
    hold off;
    
end