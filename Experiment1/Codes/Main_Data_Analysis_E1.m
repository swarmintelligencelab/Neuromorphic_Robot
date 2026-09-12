%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/06/2026
%Description : This is the main code used to perform all calculations and generate the figures for the first experiment presented in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars;
clc;
close all;

% true  = load the saved Results MAT file and only recreate figures
% false = run the complete analysis and update the MAT
useCachedResults = false;

analysisFolder = fileparts(mfilename('fullpath'));
addpath(analysisFolder);

outputFolder = fullfile(analysisFolder, 'Output_E1');
figuresFolder = fullfile(outputFolder, 'Figures');

if ~exist(figuresFolder, 'dir')
    mkdir(figuresFolder);
end

resultsFile = fullfile(outputFolder, ...
    'Result_E1.mat');

cacheFile = resultsFile;

if useCachedResults

    assert(isfile(cacheFile), ...
        'Cannot find cached result file:\n%s', cacheFile);

    load(cacheFile);

    useCachedResults = true;
    
    analysisFolder = fileparts(mfilename('fullpath'));
    outputFolder   = fullfile(analysisFolder, 'Output_E1');
    figuresFolder  = fullfile(outputFolder, 'Figures');
    resultsFile    = fullfile(outputFolder, ...
        'Result_E1.mat');
    cacheFile      = resultsFile;
    
    fprintf('Loaded cached results from:\n%s\n', cacheFile);

else

    %% 1) Load Data
    load( ...
        fullfile(analysisFolder, 'Data_E1.mat'), ...
        'Dat_contr', 'Dat_A2C', 'Dat_C2A');

    Dat_Control = Dat_contr;

    %% 2) Load Parameters
    Parameters = load_constants(Dat_Control);

    x_max = Parameters.x_max;
    y_tank = Parameters.y_tank;
    y_max = Parameters.y_max;
    fps = Parameters.fps;
    dt = Parameters.dt;
    Xmin = Parameters.Xmin;
    Xmax = Parameters.Xmax;
    Ymin = Parameters.Ymin;
    Ymax = Parameters.Ymax;
    n_trials = Parameters.n_trials;

    %% 3) Load Time series
    [XX_control, YY_control, VV_control, WW_control, ...
        AA_control, VVx_control, VVy_control] = ...
        Load_TimeSeries_Control(Parameters, Dat_Control);

    [TimeSeries_Fish_A2C, TimeSeries_Robot_A2C] = ...
        Load_TimeSeries_A_and_C(Parameters, Dat_A2C);

    [TimeSeries_Fish_C2A, TimeSeries_Robot_C2A] = ...
        Load_TimeSeries_A_and_C(Parameters, Dat_C2A);

end

set(groot, 'DefaultFigureWindowStyle', 'normal');
set(groot, ...
    'defaultAxesFontName', 'Helvetica', ...
    'defaultAxesFontSize', 12, ...
    'defaultTextFontName', 'Helvetica', ...
    'defaultLegendFontName', 'Helvetica');

pale_gray  = [0.80,0.80,0.80];
pale_black = [0.15,0.15,0.15];
pale_blue = [0.2 0.6 1];
pale_red = [1 0.4 0.4];
Font_Size = 15;
Num_Time_Bins = 2;

% Distance from the bottom
[M_YY_control] = Compute_time_Average(YY_control+y_max/2,Num_Time_Bins);
[M_YY_f_A2C] = Compute_time_Average(TimeSeries_Fish_A2C.YY_f+y_max/2,Num_Time_Bins);
[M_YY_f_C2A] = Compute_time_Average(TimeSeries_Fish_C2A.YY_f+y_max/2,Num_Time_Bins);
plot_grouped_bar(M_YY_control, M_YY_f_A2C, M_YY_f_C2A,Font_Size);
set(gcf, 'Position', [680   530   327   348]);

% Linear Speed
[M_VV_control] = Compute_time_Average(VV_control,Num_Time_Bins);
[M_VV_f_A2C] = Compute_time_Average(TimeSeries_Fish_A2C.VV_f,Num_Time_Bins);
[M_VV_f_C2A] = Compute_time_Average(TimeSeries_Fish_C2A.VV_f,Num_Time_Bins);
plot_grouped_bar(M_VV_control, M_VV_f_A2C, M_VV_f_C2A,Font_Size);
set(gcf, 'Position', [680   530   327   348]);
axis([0.5,2.5,0,30])

% Absolute Acceleration
if ~useCachedResults
    M_AA_control = Compute_time_Average( ...
        abs(AA_control), Num_Time_Bins);

    M_AA_f_A2C = Compute_time_Average( ...
        abs(TimeSeries_Fish_A2C.AA_f), Num_Time_Bins);

    M_AA_f_C2A = Compute_time_Average( ...
        abs(TimeSeries_Fish_C2A.AA_f), Num_Time_Bins);
end

plot_grouped_bar(M_AA_control, M_AA_f_A2C, M_AA_f_C2A,Font_Size);
set(gcf, 'Position', [680   530   327   348]);
axis([0.5,2.5,0,500])

% Absolute angular speed
[M_WW_control] = Compute_time_Average(abs(WW_control),Num_Time_Bins);
[M_WW_f_A2C] = Compute_time_Average(abs(TimeSeries_Fish_A2C.WW_f),Num_Time_Bins);
[M_WW_f_C2A] = Compute_time_Average(abs(TimeSeries_Fish_C2A.WW_f),Num_Time_Bins);
plot_grouped_bar(M_WW_control, M_WW_f_A2C, M_WW_f_C2A,Font_Size);
set(gcf, 'Position', [680   530   327   348]);


% Plot Robot Trajectories
Pos_R_A2C = TimeSeries_Robot_A2C.YY_r;
Pos_R_C2A = TimeSeries_Robot_C2A.YY_r;

Vel_R_A2C = TimeSeries_Robot_A2C.VV_r;
Vel_R_C2A = TimeSeries_Robot_C2A.VV_r;

Acc_R_A2C = TimeSeries_Robot_A2C.AA_r;
Acc_R_C2A = TimeSeries_Robot_C2A.AA_r;

Num_Time_Bins = 10;
[PosY_avg_A2C]  = Compute_time_Average(Pos_R_A2C+y_max/2,Num_Time_Bins);
[PosY_avg_C2A]  = Compute_time_Average(Pos_R_C2A+y_max/2,Num_Time_Bins);

[Vel_avg_A2C]   = Compute_time_Average(abs(Vel_R_A2C),Num_Time_Bins);
[Vel_avg_C2A]   = Compute_time_Average(abs(Vel_R_C2A),Num_Time_Bins);


[Accel_avg_A2C] = Compute_time_Average(abs(Acc_R_A2C),Num_Time_Bins);
[Accel_avg_C2A] = Compute_time_Average(abs(Acc_R_C2A),Num_Time_Bins);

% Position along the water column
time = 1:size(PosY_avg_A2C, 1);

mu_pos_A2C = mean(PosY_avg_A2C, 2);
mu_pos_C2A = mean(PosY_avg_C2A, 2);

sd_pos_A2C  = std(PosY_avg_A2C, 0, 2);
sd_pos_C2A  = std(PosY_avg_C2A, 0, 2);
sem_pos_A2C = sd_pos_A2C ./ sqrt(n_trials);
sem_pos_C2A = sd_pos_C2A ./ sqrt(n_trials);

figure; set(gcf, 'Position', [615 595 466*2 190*2]); hold on;
fill([time fliplr(time)], [(mu_pos_A2C+sd_pos_A2C)' fliplr((mu_pos_A2C-sd_pos_A2C)')], ...
     pale_red, 'FaceAlpha',0.20, 'EdgeColor','none');
fill([time fliplr(time)], [(mu_pos_C2A+sd_pos_C2A)' fliplr((mu_pos_C2A-sd_pos_C2A)')], ...
     pale_blue, 'FaceAlpha',0.20, 'EdgeColor','none');

errorbar(time, mu_pos_A2C, sem_pos_A2C, 'o-', 'Color', pale_red, ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 5);
errorbar(time, mu_pos_C2A, sem_pos_C2A, 'o-', 'Color', pale_blue, ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 5);

set(gca, 'FontSize',18, 'TickLabelInterpreter','latex','XColor','k','YColor','k');
set(gcf, 'Color','w');
axis([1 10 10 30]);
legend({'A-to-C','C-to-A'}, 'Location','best');
hold off;

% Linear Speed
time = 1:size(Vel_avg_A2C, 1);

mu_vel_A2C = mean(Vel_avg_A2C, 2);
mu_vel_C2A = mean(Vel_avg_C2A, 2);

sd_vel_A2C  = std(Vel_avg_A2C, 0, 2);
sd_vel_C2A  = std(Vel_avg_C2A, 0, 2);
sem_vel_A2C = sd_vel_A2C ./ sqrt(n_trials);
sem_vel_C2A = sd_vel_C2A ./ sqrt(n_trials);

figure; set(gcf, 'Position', [615 595 466*2 190*2]); hold on;
fill([time fliplr(time)], [(mu_vel_A2C+sd_vel_A2C)' fliplr((mu_vel_A2C-sd_vel_A2C)')], ...
     pale_red, 'FaceAlpha',0.20, 'EdgeColor','none');
fill([time fliplr(time)], [(mu_vel_C2A+sd_vel_C2A)' fliplr((mu_vel_C2A-sd_vel_C2A)')], ...
     pale_blue, 'FaceAlpha',0.20, 'EdgeColor','none');

errorbar(time, mu_vel_A2C, sem_vel_A2C, 'o-', 'Color', pale_red, ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 5);
errorbar(time, mu_vel_C2A, sem_vel_C2A, 'o-', 'Color', pale_blue, ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 5);

set(gca, 'FontSize',18, 'TickLabelInterpreter','latex','XColor','k','YColor','k');
set(gcf, 'Color','w');
axis([1 10 4 12]);
legend({'A-to-C','C-to-A'}, 'Location','best');
hold off;

% Linear Acceleration
time = 1:size(Accel_avg_A2C, 1);

mu_acc_A2C = mean(Accel_avg_A2C, 2);
mu_acc_C2A = mean(Accel_avg_C2A, 2);

sd_acc_A2C  = std(Accel_avg_A2C, 0, 2);
sd_acc_C2A  = std(Accel_avg_C2A, 0, 2);
sem_acc_A2C = sd_acc_A2C ./ sqrt(n_trials);
sem_acc_C2A = sd_acc_C2A ./ sqrt(n_trials);

figure; set(gcf, 'Position', [615 595 466*2 190*2]); hold on;
fill([time fliplr(time)], [(mu_acc_A2C+sd_acc_A2C)' fliplr((mu_acc_A2C-sd_acc_A2C)')], ...
     pale_red, 'FaceAlpha',0.20, 'EdgeColor','none');
fill([time fliplr(time)], [(mu_acc_C2A+sd_acc_C2A)' fliplr((mu_acc_C2A-sd_acc_C2A)')], ...
     pale_blue, 'FaceAlpha',0.20, 'EdgeColor','none');

errorbar(time, mu_acc_A2C, sem_acc_A2C, 'o-', 'Color', pale_red, ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 5);
errorbar(time, mu_acc_C2A, sem_acc_C2A, 'o-', 'Color', pale_blue, ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 5);

set(gca, 'FontSize',18, 'TickLabelInterpreter','latex','XColor','k','YColor','k');
set(gcf, 'Color','w');
axis([1 10 50 130]);
legend({'A-to-C','C-to-A'}, 'Location','best');
hold off;

% Neuron population firing
% Raster plots of neuromorphic circuit activation for all trials
Dat1 = TimeSeries_Robot_A2C.NN_raw';
Dat2 = TimeSeries_Robot_C2A.NN_raw';
[n_neurons, n_steps] = size(Dat1);
time = ((0:n_steps-1) * dt)/60;

% Shared black-to-red scale: low voltage = black, high voltage = red.
n_color_levels = 256;
black_to_red = [linspace(0,1,n_color_levels)', ...
                zeros(n_color_levels,1), ...
                zeros(n_color_levels,1)];
raw_color_max = max([Dat1(:); Dat2(:)], [], 'omitnan');
if ~isfinite(raw_color_max) || raw_color_max <= 0
    raw_color_max = 1;
end
raw_color_ticks = 0:floor(raw_color_max);

figure('Position',[680 530 327*2 348],'Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% A-to-C
ax1 = nexttile;
imagesc(ax1, time, 1:n_neurons, Dat1);
set(ax1,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex');
set(ax1,'XTickLabel',[]);
axis(ax1,[0,10,1,10]);
clim(ax1,[0 raw_color_max]);
colormap(ax1,black_to_red);
cb1 = colorbar(ax1);
set(cb1,'TickLabelInterpreter','latex','FontSize',Font_Size-2, ...
         'Ticks',raw_color_ticks,'FontName','Times New Roman');

% C-to-A
ax2 = nexttile;
imagesc(ax2, time, 1:n_neurons, Dat2);
set(ax2,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex');
axis(ax2,[0,10,1,10]);
clim(ax2,[0 raw_color_max]);
colormap(ax2,black_to_red);
cb2 = colorbar(ax2);
set(cb2,'TickLabelInterpreter','latex','FontSize',Font_Size-2, ...
         'Ticks',raw_color_ticks,'FontName','Times New Roman');

% Binary firing outputs obtained by thresholding the circuits output voltage across trials
Dat1 = TimeSeries_Robot_A2C.NN_firing';
Dat2 = TimeSeries_Robot_C2A.NN_firing';

% Shared white-to-purple scale: 0 = white, 1 = firing event.
binary_purple = [0.36 0.18 0.55];
white_to_purple = [linspace(1,binary_purple(1),n_color_levels)', ...
                   linspace(1,binary_purple(2),n_color_levels)', ...
                   linspace(1,binary_purple(3),n_color_levels)'];

figure('Position',[680 530 327*2 348],'Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% A-to-C
ax1 = nexttile;
h1 = imagesc(ax1,time,1:n_neurons,Dat1);
set(h1,'AlphaData',~isnan(Dat1));
set(ax1,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex','Color','w');
set(ax1,'XTickLabel',[]);
axis(ax1,[0,10,1,10]);
clim(ax1,[0 1]);
colormap(ax1,white_to_purple);
cb1 = colorbar(ax1);
set(cb1,'TickLabelInterpreter','latex','FontSize',Font_Size-2,...
         'Ticks',[0 1],'FontName','Times New Roman');

% C-to-A
ax2 = nexttile;
h2 = imagesc(ax2,time,1:n_neurons,Dat2);
set(h2,'AlphaData',~isnan(Dat2));
set(ax2,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex','Color','w');
axis(ax2,[0,10,1,10]);
clim(ax2,[0 1]);
colormap(ax2,white_to_purple);
cb2 = colorbar(ax2);
set(cb2,'TickLabelInterpreter','latex','FontSize',Font_Size-2,...
         'Ticks',[0 1],'FontName','Times New Roman');

Dat1 = TimeSeries_Robot_A2C.NN_firing';
Dat2 = TimeSeries_Robot_C2A.NN_firing';
[n_neurons, n_steps] = size(Dat1);
time = ((0:n_steps-1) * dt)/60;

% Neuron Firing Rate
NN_r_fire_C2A = TimeSeries_Robot_C2A.NN_firing;
NN_r_fire_A2C = TimeSeries_Robot_A2C.NN_firing;

samples_5min = round(5*60/dt);
matrices = {NN_r_fire_A2C, NN_r_fire_C2A};

N_first5min = cell(size(matrices));
N_after5min = cell(size(matrices));
freq_first5min = cell(size(matrices));
freq_after5min = cell(size(matrices));
Exposure_first5min = cell(size(matrices));
Exposure_after5min = cell(size(matrices));

for m = 1:length(matrices)

    binary_data = (matrices{m} == 1);

    % Detect firing-event onset in the FULL time series
    event_onset = [binary_data(1,:); diff(binary_data,1,1) == 1];

    % Pre
    N_first5min{m} = sum(event_onset(1:samples_5min,:),1);
    Exposure_first5min{m} = samples_5min * dt;
    freq_first5min{m} = N_first5min{m} ./ Exposure_first5min{m};

    % Post
    N_after5min{m} = sum(event_onset(samples_5min+1:end,:),1);
    Exposure_after5min{m} = size(binary_data,1)*dt - Exposure_first5min{m};
    freq_after5min{m} = N_after5min{m} ./ Exposure_after5min{m};

end

freq_neuron_A2C = [freq_first5min{1}; freq_after5min{1}];
freq_neuron_C2A = [freq_first5min{2}; freq_after5min{2}];

N_neuron_A2C = [N_first5min{1}; N_after5min{1}];
N_neuron_C2A = [N_first5min{2}; N_after5min{2}];

Exposure_neuron_A2C = [
    repmat(Exposure_first5min{1},1,size(freq_neuron_A2C,2));
    repmat(Exposure_after5min{1},1,size(freq_neuron_A2C,2))
];

Exposure_neuron_C2A = [
    repmat(Exposure_first5min{2},1,size(freq_neuron_C2A,2));
    repmat(Exposure_after5min{2},1,size(freq_neuron_C2A,2))
];

plot_grouped_2_bar(freq_neuron_A2C,freq_neuron_C2A,Font_Size)
set(gcf,'Position',[680 530 327 348]);
set(gca,'XColor','k','YColor','k');

% Linear Speed of the robot
[M_VV_r_A2C] = Compute_time_Average(TimeSeries_Robot_A2C.VV_r,2);
[M_VV_r_C2A] = Compute_time_Average(TimeSeries_Robot_C2A.VV_r,2);
plot_grouped_2_bar(M_VV_r_A2C, M_VV_r_C2A, Font_Size)
set(gcf, 'Position', [680   530   327   348]);
axis([0.5,2.5,0,15])

% Acceleration of the robot
[M_AA_r_A2C] = Compute_time_Average(abs(TimeSeries_Robot_A2C.AA_r),2);
[M_AA_r_C2A] = Compute_time_Average(abs(TimeSeries_Robot_C2A.AA_r),2);
plot_grouped_2_bar(M_AA_r_A2C, M_AA_r_C2A, Font_Size)
set(gcf, 'Position', [680   530   327   348]);
axis([0.5,2.5,0,120])

% Robot Interaction Episodes
if ~useCachedResults

    BL = 3;
    proximity_factor = 3;

    [freq_matrix_C2A, mean_dur_matrix_C2A, ...
        total_time_matrix_C2A, freq_matrix_A2C, ...
        mean_dur_matrix_A2C, total_time_matrix_A2C] = ...
        FrequencyDuration_InspectionEpisodes( ...
        TimeSeries_Fish_C2A, TimeSeries_Fish_A2C, ...
        TimeSeries_Robot_C2A, TimeSeries_Robot_A2C, ...
        dt, BL, proximity_factor);

    rng(0);
    M_chance = 10000;
    R_threshold = BL * proximity_factor;

    [chance_freq_A2C, chance_meanDur_A2C, ...
        chance_totalTime_A2C, freq_A2C, dur_A2C, total_A2C] = ...
        ComputeChanceFre_Dur_T_ShuffledGroup( ...
        TimeSeries_Fish_A2C.XX_f, TimeSeries_Fish_A2C.YY_f, ...
        TimeSeries_Robot_A2C.XX_r, TimeSeries_Robot_A2C.YY_r, ...
        dt, R_threshold, M_chance);

    [chance_freq_C2A, chance_meanDur_C2A, ...
        chance_totalTime_C2A, freq_C2A, dur_C2A, total_C2A] = ...
        ComputeChanceFre_Dur_T_ShuffledGroup( ...
        TimeSeries_Fish_C2A.XX_f, TimeSeries_Fish_C2A.YY_f, ...
        TimeSeries_Robot_C2A.XX_r, TimeSeries_Robot_C2A.YY_r, ...
        dt, R_threshold, M_chance);

    T_null = table( ...
        freq_A2C(:,1), freq_A2C(:,2), ...
        freq_C2A(:,1), freq_C2A(:,2), ...
        dur_A2C(:,1), dur_A2C(:,2), ...
        dur_C2A(:,1), dur_C2A(:,2), ...
        total_A2C(:,1), total_A2C(:,2), ...
        total_C2A(:,1), total_C2A(:,2), ...
        'VariableNames', { ...
        'Frequency_A2C_Pre','Frequency_A2C_Post', ...
        'Frequency_C2A_Pre','Frequency_C2A_Post', ...
        'MeanDuration_A2C_Pre','MeanDuration_A2C_Post', ...
        'MeanDuration_C2A_Pre','MeanDuration_C2A_Post', ...
        'TotalTime_A2C_Pre','TotalTime_A2C_Post', ...
        'TotalTime_C2A_Pre','TotalTime_C2A_Post'});

    writetable( ...
    T_null, ...
    fullfile(outputFolder, ...
    'Chance_Interaction_Null_ByCondition.xlsx'), ...
    'Sheet', 'NullDistributions');

end

% Frequency
plot_grouped_2_bar(freq_matrix_A2C, freq_matrix_C2A, Font_Size);
set(gcf, 'Position', [680 530 327 348]);
%yl = ylim; hold on;
%yline(chance_freq, '--k', 'LineWidth', 1.5);
%ylim(yl); hold off;

% Duration
plot_grouped_2_bar(mean_dur_matrix_A2C, mean_dur_matrix_C2A, Font_Size);
set(gcf, 'Position', [680 530 327 348]);
%yl = ylim; hold on;
%yline(chance_meanDur, '--k', 'LineWidth', 1.5);
%ylim(yl); hold off;

% Total time
plot_grouped_2_bar(total_time_matrix_A2C, total_time_matrix_C2A, Font_Size);
set(gcf, 'Position', [680 530 327 348]);
%yl = ylim; hold on;
%yline(chance_totalTime, '--k', 'LineWidth', 1.5);
%ylim(yl); hold off;

% Transfer Entropy Analysis
if ~useCachedResults
    % Time series
    AA_f_C2A = TimeSeries_Fish_C2A.AA_f;
    AA_f_A2C = TimeSeries_Fish_A2C.AA_f;
    
    AA_r_C2A = TimeSeries_Robot_C2A.AA_r;
    AA_r_A2C = TimeSeries_Robot_A2C.AA_r;
    
    remove_nan = 30078+1;
    AA_f_C2A(remove_nan:end,:) = [];
    AA_r_C2A(remove_nan:end,:) = [];
    AA_f_A2C(remove_nan:end,:) = [];
    AA_r_A2C(remove_nan:end,:) = [];
    AA_control_TE = AA_control(1:remove_nan-1, :);
    NN_r_fire_C2A(remove_nan:end,:) = [];
    NN_r_fire_A2C(remove_nan:end,:) = [];
    
    % Normalize signals to unit variance
    AA_f_C2A_normal = AA_f_C2A / std(AA_f_C2A(:));
    AA_f_A2C_normal = AA_f_A2C / std(AA_f_A2C(:));
    AA_r_C2A_normal = AA_r_C2A / std(AA_r_C2A(:));
    AA_r_A2C_normal = AA_r_A2C / std(AA_r_A2C(:));
    
    % Binning Time Series
    e_max = 0.06;
    nbin = 5;
    [B_f_A2C,B_r_A2C]= Bining_TimeSeries(AA_f_A2C,AA_r_A2C,e_max,nbin);
    [B_f_C2A,B_r_C2A]= Bining_TimeSeries(AA_f_C2A,AA_r_C2A,e_max,nbin);
    [B_f_control,~] = Bining_TimeSeries(AA_control_TE, AA_control_TE, e_max, nbin);
    
    % we want the symbol distrbutions to be more like uniform distributions
    % because it maximizes information content and transition diversity.
    % However, if there is noise like in our case, the central bin should be
    % higher than other because when computing the transfer entropy this bit
    % with a higher number of symbols would weight less on the TE computation
    % therefore, we would be neglecting noise and giving more weight to actual
    % data
    
    % Define parameter values
    delay_values = 0:1:10;     % Range of number of bins
    NetTE_results = zeros(1,length(delay_values));  % Preallocate
    
    for i = 1:length(delay_values)
            TAU = delay_values(i);
    
            % Compute Net TE for both conditions
            NetTE_C2A = NetTransferEntropy_main(B_f_C2A,B_r_C2A,nbin,TAU);
            NetTE_A2C = NetTransferEntropy_main(B_f_A2C, B_r_A2C, nbin, TAU);
            NetTE_results(i) = mean([mean(NetTE_C2A), mean(NetTE_A2C)]);
    end
    
    % Find best combination (maximum Net TE)
    [max_TE, idx] = max(NetTE_results(:));   
    
    
    TAU = delay_values(idx);
    
    NetTE_C2A = NetTransferEntropy_main(B_f_C2A,B_r_C2A,nbin,TAU);
    NetTE_A2C = NetTransferEntropy_main(B_f_A2C, B_r_A2C, nbin, TAU);
    
    rng(0);
    M = 10000;
    %chance_TE = ComputeChanceTE(B_f_control,B_f_C2A,B_r_C2A,B_f_A2C, B_r_A2C,nbin,TAU,M);
    %mean_chance_TE = mean(chance_TE);
    [chance_TE_A2C,TE_null_A2C] = ComputeChanceNetTE_ShuffledGroup(B_f_A2C,B_r_A2C,nbin,TAU,M);
    [chance_TE_C2A,TE_null_C2A] = ComputeChanceNetTE_ShuffledGroup(B_f_C2A,B_r_C2A,nbin,TAU,M);
end

figure;
set(gcf, 'Color', 'w');
plot(delay_values,NetTE_results) 
set(gca, 'FontSize', Font_Size, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'normal');

figure;
set(gcf, 'Position', [680 530 327*2 500], 'Color', 'w');
histogram(TE_null_A2C,50,'FaceColor',[1.0 0.4 0.4],'FaceAlpha',0.45,'EdgeColor','none','Normalization','probability'); hold on;
histogram(TE_null_C2A,50,'FaceColor',[0.2 0.6 1.0],'FaceAlpha',0.45,'EdgeColor','none','Normalization','probability');
xline(chance_TE_A2C,'--','Color',[1.0 0.4 0.4],'LineWidth',2);
xline(chance_TE_C2A,'--','Color',[0.2 0.6 1.0],'LineWidth',2);
xline(mean(NetTE_A2C),'-','Color',[1.0 0.4 0.4],'LineWidth',3);
xline(mean(NetTE_C2A),'-','Color',[0.2 0.6 1.0],'LineWidth',3);
xlabel('Net Transfer Entropy','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
set(gca,'FontSize',18,'TickLabelInterpreter','latex','LineWidth',1.2);
xlim([1.5e-3 3.7e-3]);
box off;

%figure;
%hold on;

% Plot the histogram
%histogram((chance_TE), 50, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k', 'Normalization', 'probability');

% % Plot the observed TE (mv) as a vertical red line
% set(gcf, 'Position', [680 530 327*2 500]);
% y_limits = ylim; % Get y-axis limits to scale the vertical line
% mv1 = mean(NetTE_C2A)
% mv2 = mean(NetTE_A2C)
% y_limits = [0,0.2];
% plot(([mv1 mv1]), y_limits, 'b', 'LineWidth', 2);
% plot(([mv2 mv2]), y_limits, 'r', 'LineWidth', 2);
% 
% % Formatting
% set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex', ...
%          'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.2);
% set(gcf, 'Color', 'w');
% xlim([0.7*0.001 4*0.001]);
% hold off;
% 
% [h1,p1,~,stats1] = ttest(NetTE_A2C,NetTE_C2A);
% [h2,p2,~,stats2] = ttest(NetTE_C2A,mean_chance_TE);
% [h3,p3,~,stats3] = ttest(NetTE_A2C,mean_chance_TE);
% 
% Pp = [p1,p2,p3]

A2C = NetTE_A2C(:);
C2A = NetTE_C2A(:);
Value = [C2A; A2C];
Group = [repmat("C2A", numel(C2A), 1);
         repmat("A2C", numel(A2C), 1)];  % or categorical(Group)

long_format_table_TE = table(Value, Group, ...
    'VariableNames', {'TE','Condition'});

% Save Data for R
if ~useCachedResults

save_data_excel(M_YY_control,M_YY_f_C2A,M_YY_f_A2C,fullfile(outputFolder,'Data_for_R_PositionY.xlsx'))
save_data_excel(M_VV_control,M_VV_f_C2A,M_VV_f_A2C,fullfile(outputFolder,'Data_for_R_LinearSpeed.xlsx'))
save_data_excel(M_AA_control,M_AA_f_C2A,M_AA_f_A2C,fullfile(outputFolder,'Data_for_R_LinearAcceleration.xlsx'))
save_data_excel(M_WW_control,M_WW_f_C2A,M_WW_f_A2C,fullfile(outputFolder,'Data_for_R_AngularSpeed.xlsx'))
save_data_excel_2_files(freq_matrix_C2A,freq_matrix_A2C,fullfile(outputFolder,'Data_for_R_Frequency.xlsx'))
save_data_excel_2_files(mean_dur_matrix_C2A,mean_dur_matrix_A2C,fullfile(outputFolder,'Data_for_R_Duration.xlsx'))
save_data_excel_2_files(total_time_matrix_C2A,total_time_matrix_A2C,fullfile(outputFolder,'Data_for_R_Totaltime.xlsx'))
% save_data_excel_2_files(freq_neuron_C2A,freq_neuron_A2C,fullfile(outputFolder,'Data_for_R_NuronFrequency.xlsx'))
save_data_excel_2_files(M_VV_r_C2A,M_VV_r_A2C,fullfile(outputFolder,'Data_for_R_Speed_Robot.xlsx'))
save_data_excel_2_files(M_AA_r_C2A,M_AA_r_A2C,fullfile(outputFolder,'Data_for_R_Acc_Robot.xlsx'))

% Save Neuron Firing Data for R
n_A2C = size(freq_neuron_A2C,2);
n_C2A = size(freq_neuron_C2A,2);

Subject = [(1:n_A2C)'; (1:n_A2C)'; (1:n_C2A)'; (1:n_C2A)'];
Condition = [repmat("A2C",2*n_A2C,1); repmat("C2A",2*n_C2A,1)];
Time = [ones(n_A2C,1); 2*ones(n_A2C,1); ones(n_C2A,1); 2*ones(n_C2A,1)];

Value = [freq_neuron_A2C(1,:)'; freq_neuron_A2C(2,:)'; ...
         freq_neuron_C2A(1,:)'; freq_neuron_C2A(2,:)'];

N = [N_neuron_A2C(1,:)'; N_neuron_A2C(2,:)'; ...
     N_neuron_C2A(1,:)'; N_neuron_C2A(2,:)'];

ExposureTime = [Exposure_neuron_A2C(1,:)'; Exposure_neuron_A2C(2,:)'; ...
                Exposure_neuron_C2A(1,:)'; Exposure_neuron_C2A(2,:)'];

T_neuron = table(Subject,Condition,Time,Value,N,ExposureTime);

writetable(T_neuron, fullfile(outputFolder, 'Data_for_R_NeuronFiring_Rate.xlsx'));

writetable(long_format_table_TE, fullfile(outputFolder, 'Data_for_R_TE.xlsx'));

T_TE_null = table(TE_null_A2C(:), TE_null_C2A(:), 'VariableNames', {'TE_null_A2C', 'TE_null_C2A'});

writetable(T_TE_null, fullfile(outputFolder, 'Chance_TE_Shuffled.csv'));

T_TE_summary = table(mean(NetTE_A2C), mean(NetTE_C2A), chance_TE_A2C, chance_TE_C2A, M, 'VariableNames', {'ObservedMean_A2C', 'ObservedMean_C2A', 'ChanceMean_A2C', 'ChanceMean_C2A', 'NPermutations'});

writetable(T_TE_summary, fullfile(outputFolder, 'Chance_TE_Shuffled_Summary.csv'));

size(TE_null_A2C)
size(TE_null_C2A)

T_TE_summary

end

% Interaction-metric shuffled null distributions
interactionNullData = { ...
    freq_A2C(:,1),  freq_C2A(:,1),  freq_A2C(:,2),  freq_C2A(:,2); ...
    dur_A2C(:,1),   dur_C2A(:,1),   dur_A2C(:,2),   dur_C2A(:,2); ...
    total_A2C(:,1), total_C2A(:,1), total_A2C(:,2), total_C2A(:,2)};

interactionObservedMeans = [ ...
    mean(freq_matrix_A2C(1,:),'omitnan'), ...
    mean(freq_matrix_C2A(1,:),'omitnan'), ...
    mean(freq_matrix_A2C(2,:),'omitnan'), ...
    mean(freq_matrix_C2A(2,:),'omitnan'); ...
    mean(mean_dur_matrix_A2C(1,:),'omitnan'), ...
    mean(mean_dur_matrix_C2A(1,:),'omitnan'), ...
    mean(mean_dur_matrix_A2C(2,:),'omitnan'), ...
    mean(mean_dur_matrix_C2A(2,:),'omitnan'); ...
    mean(total_time_matrix_A2C(1,:),'omitnan'), ...
    mean(total_time_matrix_C2A(1,:),'omitnan'), ...
    mean(total_time_matrix_A2C(2,:),'omitnan'), ...
    mean(total_time_matrix_C2A(2,:),'omitnan')];

plot_interaction_chance_null_grid( ...
    interactionNullData, ...
    interactionObservedMeans, ...
    Font_Size);

% Automatically save all figures
analysisFolder = fileparts(mfilename('fullpath'));
outputFolder = fullfile(analysisFolder, 'Output_E1');
figuresFolder = fullfile(outputFolder, 'Figures');

if ~isfolder(figuresFolder)
    [folderCreated, folderMessage] = mkdir(figuresFolder);
    assert(folderCreated, 'Cannot create figures folder:\n%s\n%s', figuresFolder, folderMessage);
end

figureNames = {
    '01_Fish_Distance_From_Bottom'
    '02_Fish_Linear_Speed'
    '03_Fish_Absolute_Linear_Acceleration'
    '04_Fish_Absolute_Angular_Speed'
    '05_Robot_Position_Along_Water_Column'
    '06_Robot_Linear_Speed_Time_Course'
    '07_Robot_Absolute_Acceleration_Time_Course'
    '08_Neuron_Raw_Output_Heatmap'
    '09_Neuron_Binary_Firing_Heatmap'
    '10_Neuron_Firing_Rate'
    '11_Robot_Linear_Speed_Pre_Post'
    '12_Robot_Absolute_Acceleration_Pre_Post'
    '13_Interaction_Frequency'
    '14_Interaction_Mean_Duration'
    '15_Interaction_Total_Time'
    '16_NetTE_Delay_Selection'
    '17_NetTE_Shuffled_Null_Distributions'
    '18_Frequency_A2C_Pre_ChanceNull'
    '19_Frequency_C2A_Pre_ChanceNull'
    '20_Frequency_A2C_Post_ChanceNull'close all
    '21_Frequency_C2A_Post_ChanceNull'
    
    '22_MeanDuration_A2C_Pre_ChanceNull'
    '23_MeanDuration_C2A_Pre_ChanceNull'
    '24_MeanDuration_A2C_Post_ChanceNull'
    '25_MeanDuration_C2A_Post_ChanceNull'
    
    '26_TotalTime_A2C_Pre_ChanceNull'
    '27_TotalTime_C2A_Pre_ChanceNull'
    '28_TotalTime_A2C_Post_ChanceNull'
    '29_TotalTime_C2A_Post_ChanceNull'
};

figures = findall(groot, 'Type', 'figure');
[~, order] = sort([figures.Number]);
figures = figures(order);

nExpectedFigures = numel(figureNames);

assert(numel(figures) == nExpectedFigures, ...
    'Expected %d figures, but found %d.', ...
    nExpectedFigures,numel(figures));

drawnow;

for k = 1:nExpectedFigures

    set(figures(k), ...
        'Name', figureNames{k}, ...
        'NumberTitle', 'off');

    exportgraphics( ...
        figures(k), ...
        fullfile(figuresFolder, [figureNames{k} '.png']), ...
        'Resolution', 300);

    savefig( ...
        figures(k), ...
        fullfile(figuresFolder, [figureNames{k} '.fig']));
end

fprintf('Saved %d PNG and %d FIG files to:\n%s\n', ...
    nExpectedFigures,nExpectedFigures,figuresFolder);

% Automatically save final MATLAB Workspace
if ~useCachedResults

% These graphics-handle variables are not needed for replotting.
clear ax1 ax2 cb1 cb2 t

save(resultsFile, '-v7.3');

fprintf('Final Workspace saved to:\n%s\n', resultsFile);

end
