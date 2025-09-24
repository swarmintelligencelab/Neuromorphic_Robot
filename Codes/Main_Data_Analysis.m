%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%Description : This is the main code for all calculations and visualizations of the figures in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear;
close all;


%% 1) Load Data
load('Data.mat','Dat_contr','Dat_A2C','Dat_C2A');
set(groot, 'DefaultFigureWindowStyle', 'normal');
set(groot, ...
    'defaultAxesFontName', 'Helvetica', ...
    'defaultAxesFontSize', 12, ...
    'defaultTextFontName', 'Helvetica', ...
    'defaultLegendFontName', 'Helvetica');

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
Ymax= Parameters.Ymax;
n_trials = Parameters.n_trials;


%% 3) Load Time series
[XX_control, YY_control,VV_control,WW_control,AA_control,VVx_control,VVy_control] = Load_TimeSeries_Control(Parameters,Dat_Control);

[TimeSeries_Fish_A2C,TimeSeries_Robot_A2C] = Load_TimeSeries_A_and_C(Parameters,Dat_A2C);

[TimeSeries_Fish_C2A,TimeSeries_Robot_C2A] = Load_TimeSeries_A_and_C(Parameters,Dat_C2A);

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
[M_AA_control] = Compute_time_Average(abs(AA_control),Num_Time_Bins);
[M_AA_f_A2C] = Compute_time_Average(abs(TimeSeries_Fish_A2C.AA_f),Num_Time_Bins);
[M_AA_f_C2A] = Compute_time_Average(abs(TimeSeries_Fish_C2A.AA_f),Num_Time_Bins);
plot_grouped_bar(M_AA_control, M_AA_f_A2C, M_AA_f_C2A,Font_Size);
set(gcf, 'Position', [680   530   327   348]);
axis([0.5,2.5,0,500])

% Absolute angular speed
[M_WW_control] = Compute_time_Average(abs(WW_control),Num_Time_Bins);
[M_WW_f_A2C] = Compute_time_Average(abs(TimeSeries_Fish_A2C.WW_f),Num_Time_Bins);
[M_WW_f_C2A] = Compute_time_Average(abs(TimeSeries_Fish_C2A.WW_f),Num_Time_Bins);
plot_grouped_bar(M_WW_control, M_WW_f_A2C, M_WW_f_C2A,Font_Size);
set(gcf, 'Position', [680   530   327   348]);


%==================== Plot Robot Trajectories =============================
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


%==================== Neuron population firing ============================
% Raster plots of neuromorphic circuit activation for all trials
Dat1 = TimeSeries_Robot_A2C.NN_raw';
Dat2 = TimeSeries_Robot_C2A.NN_raw';
[n_neurons, n_steps] = size(Dat1);
time = ((0:n_steps-1) * dt)/60;

figure('Position',[680 530 327*2 348],'Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% A-to-C
ax1 = nexttile;
imagesc(ax1, time, 1:n_neurons, Dat1);
set(ax1,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex');
set(ax1,'XTickLabel',[]);
axis(ax1,[0,10,1,10]);
cb1 = colorbar(ax1);
set(cb1,'TickLabelInterpreter','latex','FontSize',Font_Size-2, ...
         'Ticks',[0 1 2 3],'FontName','Times New Roman');

% C-to-A
ax2 = nexttile;
imagesc(ax2, time, 1:n_neurons, Dat2);
set(ax2,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex');
axis(ax2,[0,10,1,10]);
cb2 = colorbar(ax2);
set(cb2,'TickLabelInterpreter','latex','FontSize',Font_Size-2, ...
         'Ticks',[0 1 2 3],'FontName','Times New Roman');
colormap(gcf, cmocean('thermal'));


% Binary firing outputs obtained by thresholding the circuits output voltage across trials
Dat1 = TimeSeries_Robot_A2C.NN_firing';
Dat2 = TimeSeries_Robot_C2A.NN_firing';
figure('Position',[680 530 327*2 348],'Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% A-to-C
ax1 = nexttile;
imagesc(time,1:n_neurons,Dat1);
set(ax1,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex');
set(ax1,'XTickLabel',[]);
axis(ax1,[0,10,1,10]);
cb1 = colorbar(ax1);
set(cb1,'TickLabelInterpreter','latex','FontSize',Font_Size-2,...
         'Ticks',[0 1 2 3],'FontName','Times New Roman');

% C-to-A
ax2 = nexttile;
imagesc(time,1:n_neurons,Dat2);
set(ax2,'YDir','normal','FontSize',22,'TickLabelInterpreter','latex');
axis(ax2,[0,10,1,10]);
cb2 = colorbar(ax2);
set(cb2,'TickLabelInterpreter','latex','FontSize',Font_Size-2,...
         'Ticks',[0 1 2 3],'FontName','Times New Roman');

colormap(sky);

Dat1 = TimeSeries_Robot_A2C.NN_firing';
Dat2 = TimeSeries_Robot_C2A.NN_firing';
[n_neurons, n_steps] = size(Dat1);
time = ((0:n_steps-1) * dt)/60;


%==================== Plot Neuron Firing Rate =============================
NN_r_fire_C2A = TimeSeries_Robot_C2A.NN_firing;
NN_r_fire_A2C = TimeSeries_Robot_A2C.NN_firing;

samples_5min = round(5 * 60 / dt);
matrices = {NN_r_fire_A2C, NN_r_fire_C2A};
proportions_first5min = cell(size(matrices));
proportions_after5min = cell(size(matrices));

fs = 1/dt;

for m = 1:length(matrices)
    data = matrices{m};
    % First 5 min
    first5 = data(1:samples_5min, :);
    proportions_first5min{m} = sum(first5 == 1, 1) / size(first5, 1);
    freq_first5min{m} = proportions_first5min{m} * fs;

    % After 5 min
    after5 = data(samples_5min+1:end, :);
    proportions_after5min{m} = sum(after5 == 1, 1) / size(after5, 1);
    freq_after5min{m} = proportions_after5min{m} * fs;
end

% Firing rates of the neuromorphic circuit
freq_neuron_A2C = [freq_first5min{1};freq_after5min{1}];
freq_neuron_C2A = [freq_first5min{2};freq_after5min{2}];
plot_grouped_2_bar(freq_neuron_A2C, freq_neuron_C2A, Font_Size)
set(gcf, 'Position', [680   530   327   348]);
set(gca, 'XColor', 'k', 'YColor', 'k');

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


%==================== Robot Interaction Episodes ==========================
BL = 3;
proximity_factor = 3;
[freq_matrix_C2A, mean_dur_matrix_C2A, total_time_matrix_C2A,freq_matrix_A2C, mean_dur_matrix_A2C, total_time_matrix_A2C] = ...
FrequencyDuration_InspectionEpisodes(TimeSeries_Fish_C2A, TimeSeries_Fish_A2C,TimeSeries_Robot_C2A, TimeSeries_Robot_A2C,dt, BL, proximity_factor);

% Social interactionmetrics
rng(0);
M_chance = 4000;

XX_robot_pool = [TimeSeries_Robot_A2C.XX_r, TimeSeries_Robot_C2A.XX_r];
YY_robot_pool = [TimeSeries_Robot_A2C.YY_r, TimeSeries_Robot_C2A.YY_r];

R_threshold = BL * proximity_factor;

[chance_freq, chance_meanDur, chance_totalTime] = ComputeChanceFre_Dur_T( ...
    XX_control, YY_control, XX_robot_pool, YY_robot_pool, dt, ...
    R_threshold, M_chance); % chance level

% Frequency
plot_grouped_2_bar(freq_matrix_A2C, freq_matrix_C2A, Font_Size);
set(gcf, 'Position', [680 530 327 348]);
yl = ylim; hold on;
yline(chance_freq, '--k', 'LineWidth', 1.5);
ylim(yl); hold off;

% Duration
plot_grouped_2_bar(mean_dur_matrix_A2C, mean_dur_matrix_C2A, Font_Size);
set(gcf, 'Position', [680 530 327 348]);
yl = ylim; hold on;
yline(chance_meanDur, '--k', 'LineWidth', 1.5);
ylim(yl); hold off;

% Total time
plot_grouped_2_bar(total_time_matrix_A2C, total_time_matrix_C2A, Font_Size);
set(gcf, 'Position', [680 530 327 348]);
yl = ylim; hold on;
yline(chance_totalTime, '--k', 'LineWidth', 1.5);
ylim(yl); hold off;


%==================== Transfer Entropy Analysis ===========================
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
AA_control(remove_nan:end,:) = [];
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
[B_f_control,~]= Bining_TimeSeries(AA_control,AA_control,e_max,nbin);

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

figure;
set(gcf, 'Color', 'w');
plot(delay_values,NetTE_results) 
set(gca, 'FontSize', Font_Size, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'normal');

TAU = 2;

NetTE_C2A = NetTransferEntropy_main(B_f_C2A,B_r_C2A,nbin,TAU);
NetTE_A2C = NetTransferEntropy_main(B_f_A2C, B_r_A2C, nbin, TAU);

M = 4000;
chance_TE = ComputeChanceTE(B_f_control,B_f_C2A,B_r_C2A,B_f_A2C, B_r_A2C,nbin,TAU,M);
mean_chance_TE = mean(chance_TE);

figure;
hold on;

% Plot the histogram
histogram((chance_TE), 50, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k', 'Normalization', 'probability');

% Plot the observed TE (mv) as a vertical red line
set(gcf, 'Position', [680 530 327*2 500]);
y_limits = ylim; % Get y-axis limits to scale the vertical line
mv1 = mean(NetTE_C2A)
mv2 = mean(NetTE_A2C)
y_limits = [0,0.2];
plot(([mv1 mv1]), y_limits, 'b', 'LineWidth', 2);
plot(([mv2 mv2]), y_limits, 'r', 'LineWidth', 2);

% Formatting
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex', ...
         'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.2);
set(gcf, 'Color', 'w');
xlim([0.7*0.001 4*0.001]);
hold off;

[h1,p1,~,stats1] = ttest(NetTE_A2C,NetTE_C2A);
[h2,p2,~,stats2] = ttest(NetTE_C2A,mean_chance_TE);
[h3,p3,~,stats3] = ttest(NetTE_A2C,mean_chance_TE);

Pp = [p1,p2,p3]

A2C = NetTE_A2C(:);
C2A = NetTE_C2A(:);
Value = [C2A; A2C];
Group = [repmat("C2A", numel(C2A), 1);
         repmat("A2C", numel(A2C), 1)];  % or categorical(Group)

long_format_table_TE = table(Value, Group, ...
    'VariableNames', {'TE','Condition'});


%======================== Save Data for R =================================
save_data_excel(M_YY_control,M_YY_f_C2A,M_YY_f_A2C,'Data_for_R_PositionY.xlsx')
save_data_excel(M_VV_control,M_VV_f_C2A,M_VV_f_A2C,'Data_for_R_LinearSpeed.xlsx')
save_data_excel(M_AA_control,M_AA_f_C2A,M_AA_f_A2C,'Data_for_R_LinearAcceleration.xlsx')
save_data_excel(M_WW_control,M_WW_f_C2A,M_WW_f_A2C,'Data_for_R_AngularSpeed.xlsx')
save_data_excel_2_files(freq_matrix_C2A,freq_matrix_A2C,'Data_for_R_Frequency.xlsx')
save_data_excel_2_files(mean_dur_matrix_C2A,mean_dur_matrix_A2C,'Data_for_R_Duration.xlsx')
save_data_excel_2_files(total_time_matrix_C2A,total_time_matrix_A2C,'Data_for_R_Totaltime.xlsx')
save_data_excel_2_files(freq_neuron_C2A,freq_neuron_A2C,'Data_for_R_NuronFrequency.xlsx')
save_data_excel_2_files(M_VV_r_C2A,M_VV_r_A2C,'Data_for_R_Speed_Robot.xlsx')
save_data_excel_2_files(M_AA_r_C2A,M_AA_r_A2C,'Data_for_R_Acc_Robot.xlsx')
writetable(long_format_table_TE, 'Data_for_R_TE.xlsx');