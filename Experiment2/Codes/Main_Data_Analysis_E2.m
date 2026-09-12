%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%Description : This is the main code used to perform all calculations and generate the figures for the second experiment presented in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars;
clc;
close all;

% true  = load the saved Results MAT file and only recreate figures
% false = run the complete analysis and update the MAT
useCachedResults = false;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

dataFile = fullfile(scriptDir, 'Data_E2.mat');
outputFolder = fullfile(scriptDir, 'Output_E2');
cacheFile = fullfile(outputFolder, 'Result_E2.mat');
showFigures = true;

analysisOptions = struct();
analysisOptions.ShowSignificance = true;
analysisOptions.StatsFolder = fullfile( ...
scriptDir, 'Statistics_new', 'Stats_results');

figureDir = fullfile(outputFolder, 'Figures');
if ~isfolder(outputFolder), mkdir(outputFolder); end
if ~isfolder(figureDir), mkdir(figureDir); end

if useCachedResults
    assert(isfile(cacheFile), ...
        ['Cached results file not found:\n%s\n' ...
         'Run once with useCachedResults = false, or save the current workspace first.'], ...
        cacheFile);
    cached = load(cacheFile, 'Results');
    assert(isfield(cached, 'Results'), ...
        'The cache file does not contain the Results structure: %s', cacheFile);
    Results = cached.Results;

    % Restore the variables required by the shared plotting sections.
    Parameters = Results.Parameters;
    Metadata = Results.Metadata;
    conditionDefs = Results.ConditionDefinitions;
    TrialIDs = Results.TrialIDs;
    TimeSeries_Fish = Results.TimeSeries_Fish;
    TimeSeries_Robot = Results.TimeSeries_Robot;
    nConditions = numel(conditionDefs);

    M_YY_f = Results.Fish.PositionY;
    M_VV_f = Results.Fish.LinearSpeed;
    M_AA_f = Results.Fish.LinearAcceleration;
    M_WW_f = Results.Fish.AngularSpeed;

    M_VV_r = Results.Robot.LinearSpeed;
    M_AA_r = Results.Robot.LinearAcceleration;
    RobotPosY_10Bin = Results.Robot.PositionY_10Bin;
    RobotSpeed_10Bin = Results.Robot.LinearSpeed_10Bin;
    RobotAcc_10Bin = Results.Robot.LinearAcceleration_10Bin;

    firing_rate = Results.Neuron.FiringRate;
    firing_count = Results.Neuron.EventCount;
    firing_exposure = Results.Neuron.ExposureTime;

    freq_matrix = Results.Interaction.Frequency;
    mean_dur_matrix = Results.Interaction.MeanDuration;
    total_time_matrix = Results.Interaction.TotalTime;

    fprintf('Loaded cached results:\n%s\n', cacheFile);
else
    rng(1);


%% 1) Load Data
S = load(dataFile);
conditionDefs = build_condition_definitions();
assert_required_data(S, conditionDefs);

if isfield(S, 'Metadata')
Metadata = S.Metadata;
else
Metadata = struct();
end

Parameters = load_constants_6conditions(S, Metadata, conditionDefs);

fprintf('Loaded %s\n', dataFile);
fprintf('FPS %.4g, dt %.4g, frames %d\n', Parameters.fps, Parameters.dt, Parameters.T_total);


%% 2) Load Time Series
nConditions = numel(conditionDefs);
TimeSeries_Fish = cell(1, nConditions);
TimeSeries_Robot = cell(1, nConditions);
TrialIDs = cell(1, nConditions);

for c = 1:nConditions
    datName = ['Dat_' conditionDefs(c).name];
Dat = S.(datName);
TrialIDs{c} = [Dat.TrialID];
[TimeSeries_Fish{c}, TimeSeries_Robot{c}] = load_time_series_condition(Parameters, Dat);
fprintf('%s: %d trials\n', conditionDefs(c).display, numel(Dat));
end


Results = struct();
Results.Parameters = Parameters;
Results.Metadata = Metadata;
Results.ConditionDefinitions = conditionDefs;
Results.TrialIDs = TrialIDs;
Results.TimeSeries_Fish = TimeSeries_Fish;
Results.TimeSeries_Robot = TimeSeries_Robot;
Results.SwitchWindowAvailability = summarize_switch_window_availability( ...
TimeSeries_Robot, TrialIDs, conditionDefs, Parameters);
writetable(Results.SwitchWindowAvailability, ...
fullfile(outputFolder, 'SwitchWindowAvailability.xlsx'));


%% 3) Fish behavior
numTimeBins = 2;
M_YY_f = cell(1, nConditions);
M_VV_f = cell(1, nConditions);
M_AA_f = cell(1, nConditions);
M_WW_f = cell(1, nConditions);

for c = 1:nConditions
    Fish = TimeSeries_Fish{c};
M_YY_f{c} = compute_switch_phase_average(Fish.YY_f + Parameters.y_max/2, ...
RobotPhase(TimeSeries_Robot{c}), Parameters);
M_VV_f{c} = compute_switch_phase_average(Fish.VV_f, ...
RobotPhase(TimeSeries_Robot{c}), Parameters);
M_AA_f{c} = compute_switch_phase_average(abs(Fish.AA_f), ...
RobotPhase(TimeSeries_Robot{c}), Parameters);
M_WW_f{c} = compute_switch_phase_average(abs(Fish.WW_f), ...
RobotPhase(TimeSeries_Robot{c}), Parameters);
end

Results.Fish.PositionY = M_YY_f;
Results.Fish.LinearSpeed = M_VV_f;
Results.Fish.LinearAcceleration = M_AA_f;
Results.Fish.AngularSpeed = M_WW_f;


%% 4) Robot behavior
M_VV_r = cell(1,nConditions);
M_AA_r = cell(1,nConditions);

for c = 1:nConditions
    Robot = TimeSeries_Robot{c};

    M_VV_r{c} = compute_switch_phase_average( ...
        Robot.VV_r, RobotPhase(Robot), Parameters);

    M_AA_r{c} = compute_switch_phase_average( ...
        abs(Robot.AA_r), RobotPhase(Robot), Parameters);
end

Results.Robot.LinearSpeed = M_VV_r;
Results.Robot.LinearAcceleration = M_AA_r;


%% Robot time-resolved trajectories
numRobotTimeBins = 10;

RobotPosY_10Bin = cell(1,nConditions);
RobotSpeed_10Bin = cell(1,nConditions);
RobotAcc_10Bin = cell(1,nConditions);

for c = 1:nConditions

    Robot = TimeSeries_Robot{c};

    RobotPosY_10Bin{c} = compute_switch_timebin_average( ...
        Robot.YY_r + Parameters.y_max/2, ...
        RobotPhase(Robot), ...
        Parameters, ...
        numRobotTimeBins);

    RobotSpeed_10Bin{c} = compute_switch_timebin_average( ...
        abs(Robot.VV_r), ...
        RobotPhase(Robot), ...
        Parameters, ...
        numRobotTimeBins);

    RobotAcc_10Bin{c} = compute_switch_timebin_average( ...
        abs(Robot.AA_r), ...
        RobotPhase(Robot), ...
        Parameters, ...
        numRobotTimeBins);
end

Results.Robot.PositionY_10Bin = RobotPosY_10Bin;
Results.Robot.LinearSpeed_10Bin = RobotSpeed_10Bin;
Results.Robot.LinearAcceleration_10Bin = RobotAcc_10Bin;
end


%% 5) Common plotting settings
paleColors = cat(1, conditionDefs.color);
groupedPlotOrder = grouped_condition_plot_order(conditionDefs);
groupedConditionDefs = conditionDefs(groupedPlotOrder);
groupedColors = paleColors(groupedPlotOrder, :);
fontSize = 15;
visibleState = ternary(showFigures, 'on', 'off');
phaseLabels = {'C/pre-switch 300 s', 'A/post-switch 300 s'};
tenBinLabels = compose('Bin %d', 1:10);
statsFolder = get_option(analysisOptions, 'StatsFolder', ...
fullfile(scriptDir, 'Statistics_new', 'Results1'));
showSignificance = get_option(analysisOptions, 'ShowSignificance', true);
if showSignificance && ~isfolder(statsFolder)
warning('Statistics folder not found; significance annotations disabled: %s', statsFolder);
showSignificance = false;
end


%% Plot robot 10-bin trajectories
plot_timebin_trajectories( ...
    RobotPosY_10Bin, ...
    conditionDefs, ...
    paleColors, ...
    'Robot position from bottom (cm)', ...
    'Robot Position Y', ...
    visibleState, ...
    figureDir, ...
    'Robot_PositionY_10Bin');

plot_timebin_trajectories( ...
    RobotSpeed_10Bin, ...
    conditionDefs, ...
    paleColors, ...
    'Robot linear speed (cm/s)', ...
    'Robot Linear Speed', ...
    visibleState, ...
    figureDir, ...
    'Robot_LinearSpeed_10Bin');

plot_timebin_trajectories( ...
    RobotAcc_10Bin, ...
    conditionDefs, ...
    paleColors, ...
    'Robot absolute acceleration (cm/s^2)', ...
    'Robot Absolute Acceleration', ...
    visibleState, ...
    figureDir, ...
    'Robot_Acceleration_10Bin');


%% 6) Neuron firing
if ~useCachedResults
firing_rate = cell(1,nConditions);
firing_count = cell(1,nConditions);
firing_exposure = cell(1,nConditions);

for c = 1:nConditions

    [firing_rate{c}, ...
     firing_count{c}, ...
     firing_exposure{c}] = ...
        compute_firing_event_rate( ...
            TimeSeries_Robot{c},Parameters);

end

Results.Neuron.FiringRate = firing_rate;
Results.Neuron.EventCount = firing_count;
Results.Neuron.ExposureTime = firing_exposure;
end


%% Neuron raster visualization
Results.Neuron.RasterAlgorithms = ...
    plot_neuron_raster_algorithms( ...
        TimeSeries_Robot, ...
        conditionDefs, ...
        Parameters, ...
        visibleState, ...
        figureDir);


%% Neuron firing-frequency barplot
plot_prepost_object_metric( ...
    firing_rate, ...
    'Firing frequency (Hz)', ...
    'Neuron firing frequency', ...
    visibleState, ...
    figureDir, ...
    'Neuron_FiringFrequency_PrePost');


%% Robot acceleration barplot
plot_prepost_object_metric( ...
    M_AA_r, ...
    'Absolute linear acceleration (cm/s^2)', ...
    'Robot absolute linear acceleration', ...
    visibleState, ...
    figureDir, ...
    'Robot_AbsoluteLinearAcceleration_PrePost');

if false

figure('Color','w','Position',[50 100 1600 430]);

tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

% Variables and labels
DATA = {M_YY_f, M_VV_f, M_AA_f, M_WW_f};

yLabels = {'$M_{YY}$', ...
'$M_{VV}$', ...
'$M_{AA}$', ...
'$M_{WW}$'};

panelLabels = {'(A)','(B)','(C)','(D)'};

objectNames = {'Rod','Rectangle','Fish model', ...
'Rod','Rectangle','Fish model'};

% Index order
idxIn  = [3 2 1];   % Rod, Rectangle, Fish model
idxOut = [6 5 4];   % Rod, Rectangle, Fish model

% x centers for the 6 groups:
% Inside = 1,2,3
% Outside = 5,6,7
xCenters = [1 2 3 5 6 7];

% Colors
Cpre  = [0.72 0.62 0.85];   % light purple
Cpost = [0.36 0.18 0.55];   % dark purple

% Loop over the four variables
for p = 1:4

ax = nexttile;
hold(ax,'on');

M = DATA{p};

% Build Y and SEM arrays for 6 groups x 2 bars
% columns: [Pre Post]
% rows 1:3 = Inside; rows 4:6 = Outside
Y   = zeros(6,2);
SEM = zeros(6,2);

% Inside
for j = 1:3
    D = M{idxIn(j)};
    Y(j,:) = mean(D,2,'omitnan')';

    n = sum(~isnan(D),2);
    SEM(j,:) = (std(D,0,2,'omitnan')./sqrt(n))';
end

% Outside
for j = 1:3
    D = M{idxOut(j)};
    Y(3+j,:) = mean(D,2,'omitnan')';

    n = sum(~isnan(D),2);
    SEM(3+j,:) = (std(D,0,2,'omitnan')./sqrt(n))';
end

% Draw grouped bars at custom x positions
b = bar(ax,xCenters,Y,'grouped','BarWidth',0.75);

b(1).FaceColor = Cpre;
b(1).EdgeColor = 'k';
b(1).LineWidth = 0.8;

b(2).FaceColor = Cpost;
b(2).EdgeColor = 'k';
b(2).LineWidth = 0.8;

% Error bars
for k = 1:2
    x = b(k).XEndPoints;
    errorbar(ax,x,Y(:,k),SEM(:,k), ...
        'k','LineStyle','none', ...
        'LineWidth',1.1,'CapSize',5);
end

% Vertical dashed separator between Inside and Outside
xline(ax,4,'--','Color',[0.6 0.6 0.6],'LineWidth',1.2);

% Axes formatting
xlim(ax,[0.3 7.7]);

xticks(ax,xCenters);
xticklabels(ax,objectNames);
xtickangle(ax,30);

ylabel(ax,yLabels{p}, ...
    'Interpreter','latex', ...
    'FontSize',15);

title(ax,panelLabels{p}, ...
    'FontWeight','bold', ...
    'FontSize',13);

set(ax, ...
    'FontSize',11, ...
    'LineWidth',1.0, ...
    'TickDir','out');

box(ax,'off');

% Add "Inside tank" and "Outside tank" labels
yl = ylim(ax);
yr = yl(2)-yl(1);

ylim(ax,[yl(1)-0.14*yr yl(2)]);
yl = ylim(ax);

yText = yl(1) + 0.03*(yl(2)-yl(1));

text(ax,2,yText,'Inside tank', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',10);

text(ax,6,yText,'Outside tank', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',10);

% Legend only in first panel
if p == 1
    legend(ax,{'Pre-Switch','Post-Switch'}, ...
        'Location','best', ...
        'Box','off', ...
        'FontSize',10);
end

end

end % legacy combined fish-behavior panel


%% Individual fish-behavior figures
plot_prepost_object_metric(M_YY_f, ...
    'Distance from bottom (cm)', 'Fish vertical position', ...
    visibleState, figureDir, 'Fish_PositionY_PrePost');

plot_prepost_object_metric(M_VV_f, ...
    'Linear speed (cm/s)', 'Fish linear speed', ...
    visibleState, figureDir, 'Fish_LinearSpeed_PrePost');

plot_prepost_object_metric(M_AA_f, ...
    'Absolute linear acceleration (cm/s^2)', ...
    'Fish absolute linear acceleration', ...
    visibleState, figureDir, 'Fish_AbsoluteLinearAcceleration_PrePost');

plot_prepost_object_metric(M_WW_f, ...
    'Absolute angular speed (rad/s)', 'Fish absolute angular speed', ...
    visibleState, figureDir, 'Fish_AbsoluteAngularSpeed_PrePost');


%% 7) Robot/object interaction episodes
if ~useCachedResults

BL = 3;
proximityFactor = 2;
R_threshold = BL * proximityFactor;
M_chance = 10000;


% =========================================================================
% Pool definitions for chance analysis
%
% Pool 1: Inside  = Fish + Rectangle + Rod
% Pool 2: Outside = Fish + Rectangle + Rod
% =========================================================================
poolIdx = {1:3, 4:6};
poolNames = {'In','Out'};
nPools = 2;


freq_matrix = cell(1,nConditions);
mean_dur_matrix = cell(1,nConditions);
total_time_matrix = cell(1,nConditions);

event_count_matrix = cell(1,nConditions);
observation_time_matrix = cell(1,nConditions);

% Observed interactions
for c = 1:nConditions

    [freq_matrix{c}, ...
     mean_dur_matrix{c}, ...
     total_time_matrix{c}, ...
     event_count_matrix{c}, ...
     observation_time_matrix{c}] = ...
        compute_interaction_episodes( ...
            TimeSeries_Fish{c}, ...
            TimeSeries_Robot{c}, ...
            Parameters, ...
            R_threshold);

end

% Location-pooled shuffled null distributions
chance_freq_pool = cell(1,nPools);
chance_dur_pool = cell(1,nPools);
chance_total_pool = cell(1,nPools);

null_freq_pool = cell(1,nPools);
null_dur_pool = cell(1,nPools);
null_total_pool = cell(1,nPools);

for g = 1:nPools

    idx = poolIdx{g};

    [FishPool,RobotPool] = pool_conditions_for_interaction( ...
        TimeSeries_Fish, ...
        TimeSeries_Robot, ...
        idx);

    [chance_freq_pool{g}, ...
     chance_dur_pool{g}, ...
     chance_total_pool{g}, ...
     null_freq_pool{g}, ...
     null_dur_pool{g}, ...
     null_total_pool{g}] = ...
        compute_chance_interactions_shuffled_group( ...
            FishPool, ...
            RobotPool, ...
            Parameters, ...
            R_threshold, ...
            M_chance);

end

% Save to Results
Results.Interaction.Frequency = freq_matrix;
Results.Interaction.MeanDuration = mean_dur_matrix;
Results.Interaction.TotalTime = total_time_matrix;
Results.Interaction.EventCount = event_count_matrix;
Results.Interaction.ObservationTime = observation_time_matrix;

Results.Interaction.PoolNames = poolNames;

Results.Interaction.ChanceFrequencyByLocation = chance_freq_pool;
Results.Interaction.ChanceDurationByLocation = chance_dur_pool;
Results.Interaction.ChanceTotalTimeByLocation = chance_total_pool;

Results.Interaction.NullFrequencyByLocation = null_freq_pool;
Results.Interaction.NullDurationByLocation = null_dur_pool;
Results.Interaction.NullTotalTimeByLocation = null_total_pool;


%% 8) Transfer Entropy Analysis
e_max = 0.06*10;   % intentional
nbin = 5;

tePhaseNames = {'Pre','Post'};
nPhases = 2;

% Rows = Location: In / Out
% Columns = Phase: Pre / Post
B_f_pool = cell(nPools,nPhases);
B_r_pool = cell(nPools,nPhases);

% Metadata for exporting individual-trial TE
teConditionPool = cell(nPools,1);
teTrialIDPool = cell(nPools,1);


% =========================================================================
% Pool Fish / Rectangle / Rod within each Location,
% but preserve Pre and Post separately.
% =========================================================================
for g = 1:nPools

    idxConditions = poolIdx{g};

    AA_f_phase = cell(1,nPhases);
    AA_r_phase = cell(1,nPhases);

    teConditionPool{g} = strings(0,1);
    teTrialIDPool{g} = [];

    for c = idxConditions

        Fish = TimeSeries_Fish{c};
        Robot = TimeSeries_Robot{c};

        nTrials = size(Fish.AA_f,2);

        for tr = 1:nTrials

            phaseIdx = switch_phase_indices( ...
                Robot.SwitchIdx(tr), ...
                Robot.NumFrames(tr), ...
                Parameters);

            % For paired Pre/Post TE, require both complete windows
            if isempty(phaseIdx{1}) || isempty(phaseIdx{2})
                warning( ...
                    'Skipping TE: %s, trial %g lacks a complete Pre/Post window.', ...
                    conditionDefs(c).display, ...
                    TrialIDs{c}(tr));
                continue;
            end

            trialFish = cell(1,nPhases);
            trialRobot = cell(1,nPhases);
            isValid = true;

            for phase = 1:nPhases

                trialFish{phase} = ...
                    Fish.AA_f(phaseIdx{phase},tr);

                trialRobot{phase} = ...
                    Robot.AA_r(phaseIdx{phase},tr);

                trialFish{phase} = trialFish{phase}(:);
                trialRobot{phase} = trialRobot{phase}(:);

                if any(~isfinite(trialFish{phase})) || ...
                   any(~isfinite(trialRobot{phase}))
                    isValid = false;
                end
            end

            % Do not let NaN become state zero during binning
            if ~isValid
                warning( ...
                    'Skipping TE: %s, trial %g contains non-finite acceleration.', ...
                    conditionDefs(c).display, ...
                    TrialIDs{c}(tr));
                continue;
            end

            for phase = 1:nPhases
                AA_f_phase{phase}(:,end+1) = trialFish{phase};
                AA_r_phase{phase}(:,end+1) = trialRobot{phase};
            end

            teConditionPool{g}(end+1,1) = ...
                string(conditionDefs(c).display);

            teTrialIDPool{g}(end+1,1) = ...
                TrialIDs{c}(tr);
        end
    end

    nKept = size(AA_f_phase{1},2);

    if nKept < 2
        error( ...
            'Location %s has fewer than two complete TE trials.', ...
            poolNames{g});
    end

    assert( ...
        size(AA_f_phase{2},2) == nKept, ...
        'Pre/Post TE trial counts differ.');

    % Bin Pre and Post together so both phases use exactly the same
    % epsilon/bin thresholds within this Location.
    AA_f_all = [AA_f_phase{1}, AA_f_phase{2}];
    AA_r_all = [AA_r_phase{1}, AA_r_phase{2}];

    [B_f_all,B_r_all] = Bining_TimeSeries( ...
        AA_f_all, ...
        AA_r_all, ...
        e_max, ...
        nbin);

    B_f_pool{g,1} = B_f_all(:,1:nKept);
    B_f_pool{g,2} = B_f_all(:,nKept+1:end);

    B_r_pool{g,1} = B_r_all(:,1:nKept);
    B_r_pool{g,2} = B_r_all(:,nKept+1:end);
end


% =========================================================================
% One common TAU across In/Out AND Pre/Post
% =========================================================================
delay_values = 0:10;
NetTE_results = nan(size(delay_values));

for i = 1:numel(delay_values)

    tauNow = delay_values(i);
    TE_all = [];

    for g = 1:nPools
        for phase = 1:nPhases

            TE_temp = NetTransferEntropy_main( ...
                B_f_pool{g,phase}, ...
                B_r_pool{g,phase}, ...
                nbin, ...
                tauNow);

            TE_all = [TE_all; TE_temp(:)];
        end
    end

    NetTE_results(i) = mean(TE_all,'omitnan');
end

[max_TE,tauIndex] = max(NetTE_results);

TAU = delay_values(tauIndex);

fprintf('Optimal common TAU = %d\n',TAU);
fprintf('Maximum mean TE = %.6f\n',max_TE);


% =========================================================================
% Observed TE: Location x Phase
% =========================================================================
NetTE_pool = cell(nPools,nPhases);

for g = 1:nPools
    for phase = 1:nPhases

        NetTE_pool{g,phase} = NetTransferEntropy_main( ...
            B_f_pool{g,phase}, ...
            B_r_pool{g,phase}, ...
            nbin, ...
            TAU);
    end
end


% =========================================================================
% Shuffled TE null: Location x Phase
% =========================================================================
M_TE = 10000;

[chance_TE_pool,TE_null_pool,TAU_null] = ...
    ComputeChanceNetTE_ShuffledGroup_MultiCondition( ...
        B_f_pool, ...
        B_r_pool, ...
        nbin, ...
        delay_values, ...
        M_TE);


% =========================================================================
% Results
% =========================================================================
Results.TE.PoolNames = poolNames;
Results.TE.PhaseNames = tePhaseNames;

Results.TE.NetTEByLocationPhase = NetTE_pool;
Results.TE.ChanceByLocationPhase = chance_TE_pool;
Results.TE.NullByLocationPhase = TE_null_pool;

Results.TE.TAU = TAU;
Results.TE.DelayValues = delay_values;
Results.TE.DelaySweep = NetTE_results;
Results.TE.NullTAU = TAU_null;

Results.TE.ConditionByLocation = teConditionPool;
Results.TE.TrialIDByLocation = teTrialIDPool;


%% Save Data for R
conditionNames = {'In-FishModel', 'In-Rectangle', 'In-Rod', ...
                  'Out-FishModel', 'Out-Rectangle', 'Out-Rod'};


%% Save pooled interaction null
Location = strings(0,1);
Permutation = [];
Time = [];

Frequency = [];
MeanDuration = [];
TotalTime = [];

for g = 1:nPools

    nPerm = size(null_freq_pool{g},1);

    for phase = 1:2

        Location = [ ...
            Location; ...
            repmat(string(poolNames{g}),nPerm,1) ...
        ];

        Permutation = [ ...
            Permutation; ...
            (1:nPerm)' ...
        ];

        Time = [ ...
            Time; ...
            repmat(phase,nPerm,1) ...
        ];

        Frequency = [ ...
            Frequency; ...
            null_freq_pool{g}(:,phase) ...
        ];

        MeanDuration = [ ...
            MeanDuration; ...
            null_dur_pool{g}(:,phase) ...
        ];

        TotalTime = [ ...
            TotalTime; ...
            null_total_pool{g}(:,phase) ...
        ];

    end
end


T_null = table( ...
    Location, ...
    Permutation, ...
    Time, ...
    Frequency, ...
    MeanDuration, ...
    TotalTime);


writetable( ...
    T_null, ...
    fullfile( ...
        outputFolder, ...
        'Chance_Interaction_Null_ByLocation.xlsx'), ...
    'Sheet', ...
    'NullDistributions');


%% Save observed metrics
metricData = {M_YY_f, M_VV_f, M_AA_f, M_WW_f, ...
              freq_matrix, mean_dur_matrix, total_time_matrix, ...
              M_VV_r, M_AA_r};

fileNames = {'Data_for_R_PositionY.xlsx', ...
             'Data_for_R_LinearSpeed.xlsx', ...
             'Data_for_R_LinearAcceleration.xlsx', ...
             'Data_for_R_AngularSpeed.xlsx', ...
             'Data_for_R_Frequency.xlsx', ...
             'Data_for_R_Duration.xlsx', ...
             'Data_for_R_Totaltime.xlsx', ...
             'Data_for_R_Speed_Robot.xlsx', ...
             'Data_for_R_Acc_Robot.xlsx'};


for m = 1:length(metricData)

M = metricData{m};

Condition = {};
TrialID = [];
Time = [];
Value = [];

% Only used for Frequency
EventCount = [];
TotalObservationTime = [];

for c = 1:nConditions

    D = M{c};
    nTrials = size(D,2);

    for tr = 1:nTrials

        % =================================================================
        % Pre-switch
        % =================================================================
        Condition{end+1,1} = conditionNames{c};
        TrialID(end+1,1) = TrialIDs{c}(tr);
        Time(end+1,1) = 1;
        Value(end+1,1) = D(1,tr);

        % m == 5 corresponds to Frequency
        if m == 5
            EventCount(end+1,1) = ...
                event_count_matrix{c}(1,tr);

            TotalObservationTime(end+1,1) = ...
                observation_time_matrix{c}(1,tr);
        end

        % =================================================================
        % Post-switch
        % =================================================================
        Condition{end+1,1} = conditionNames{c};
        TrialID(end+1,1) = TrialIDs{c}(tr);
        Time(end+1,1) = 2;
        Value(end+1,1) = D(2,tr);

        if m == 5
            EventCount(end+1,1) = ...
                event_count_matrix{c}(2,tr);

            TotalObservationTime(end+1,1) = ...
                observation_time_matrix{c}(2,tr);
        end

    end
end

% Build table
if m == 5

    % Frequency gets two extra columns
    T = table( ...
        Condition, ...
        TrialID, ...
        Time, ...
        Value, ...
        EventCount, ...
        TotalObservationTime);

else

    % Other metrics unchanged
    T = table( ...
        Condition, ...
        TrialID, ...
        Time, ...
        Value);
end


writetable(T, fullfile(outputFolder, fileNames{m}));

end


%% Save neuron firing data for R
Condition = {};
TrialID = [];
Time = [];
Value = [];
N = [];
ExposureTime = [];

for c = 1:nConditions

    nTrials = size(firing_rate{c},2);

    for tr = 1:nTrials
        for phase = 1:2

            Condition{end+1,1} = conditionNames{c};
            TrialID(end+1,1) = TrialIDs{c}(tr);
            Time(end+1,1) = phase;

            Value(end+1,1) = firing_rate{c}(phase,tr);
            N(end+1,1) = firing_count{c}(phase,tr);
            ExposureTime(end+1,1) = firing_exposure{c}(phase,tr);
        end
    end
end

T_neuron = table(Condition,TrialID,Time,Value,N,ExposureTime);

writetable(T_neuron, ...
    fullfile(outputFolder,'Data_for_R_NeuronFiring_Rate.xlsx'));


%% Save observed TE for R
Location = strings(0,1);
Condition = strings(0,1);
TrialID = [];
Time = [];
TE = [];

for g = 1:nPools

    nTrials = numel(teTrialIDPool{g});

    for phase = 1:nPhases

        x = NetTE_pool{g,phase}(:);

        assert( ...
            numel(x) == nTrials, ...
            'TE value count does not match metadata.');

        Location = [ ...
            Location; ...
            repmat(string(poolNames{g}),nTrials,1)];

        Condition = [ ...
            Condition; ...
            teConditionPool{g}];

        TrialID = [ ...
            TrialID; ...
            teTrialIDPool{g}];

        Time = [ ...
            Time; ...
            repmat(phase,nTrials,1)];

        TE = [TE; x];
    end
end

T_TE = table( ...
    Location, ...
    Condition, ...
    TrialID, ...
    Time, ...
    TE);

writetable( ...
    T_TE, ...
    fullfile(outputFolder,'Data_for_R_TE.xlsx'));


%% Save TE null for R
Location = strings(0,1);
Time = [];
Permutation = [];
TE_null_value = [];
TAU_null_value = [];

for g = 1:nPools
    for phase = 1:nPhases

        nPerm = numel(TE_null_pool{g,phase});

        Location = [ ...
            Location; ...
            repmat(string(poolNames{g}),nPerm,1)];

        Time = [ ...
            Time; ...
            repmat(phase,nPerm,1)];

        Permutation = [ ...
            Permutation; ...
            (1:nPerm)'];

        TE_null_value = [ ...
            TE_null_value; ...
            TE_null_pool{g,phase}(:)];

        TAU_null_value = [ ...
            TAU_null_value; ...
            TAU_null(:)];
    end
end

T_TE_null = table( ...
    Location, ...
    Time, ...
    Permutation, ...
    TE_null_value, ...
    TAU_null_value, ...
    'VariableNames', ...
    {'Location','Time','Permutation','TE_null','TAU_null'});

writetable( ...
    T_TE_null, ...
    fullfile(outputFolder,'Chance_TE_Shuffled.xlsx'));


%% Save TE summary for R
Location = strings(nPools*nPhases,1);
Time = nan(nPools*nPhases,1);
ObservedMean = nan(nPools*nPhases,1);
ChanceMean = nan(nPools*nPhases,1);
TAU_summary = repmat(TAU,nPools*nPhases,1);
NPermutations = repmat(M_TE,nPools*nPhases,1);

row = 0;

for g = 1:nPools
    for phase = 1:nPhases

        row = row + 1;

        Location(row) = string(poolNames{g});
        Time(row) = phase;

        ObservedMean(row) = ...
            mean(NetTE_pool{g,phase},'omitnan');

        ChanceMean(row) = ...
            mean(TE_null_pool{g,phase},'omitnan');
    end
end

T_TE_summary = table( ...
    Location, ...
    Time, ...
    ObservedMean, ...
    ChanceMean, ...
    TAU_summary, ...
    NPermutations);

writetable( ...
    T_TE_summary, ...
    fullfile(outputFolder,'Chance_TE_Shuffled_Summary.xlsx'));
end

% Figure: interaction metrics
if false
figure('Color','w','Position',[100 100 1250 420]);

tiledlayout(1,3, ...
'TileSpacing','compact', ...
'Padding','compact');

% Data
DATA = {freq_matrix, ...
mean_dur_matrix, ...
total_time_matrix};

yLabels = {'Frequency', ...
'Mean duration', ...
'Total interaction time'};

panelLabels = {'(A)','(B)','(C)'};

% Organization
idxIn  = [3 2 1]; % Rod, Rectangle, Fish model
idxOut = [6 5 4];

objectNames = {'Rod','Rectangle','Fish model', ...
'Rod','Rectangle','Fish model'};

% Gap between Inside and Outside
xCenters = [1 2 3 5 6 7];

% Colors
Cpre  = [0.72 0.62 0.85];   % light purple
Cpost = [0.36 0.18 0.55];   % dark purple

% Loop over metrics
for p = 1:3

ax = nexttile;
hold(ax,'on');

M = DATA{p};

% Rows = six experimental conditions
% Columns = Pre / Post
Y   = zeros(6,2);
SEM = zeros(6,2);

idx = [idxIn idxOut];

% Extract data and compute mean +/- SEM
for j = 1:6

    D = M{idx(j)};

    % =====================================================================
    % If D is 2 x N:
    % row 1 = Pre
    % row 2 = Post
    % =====================================================================
    
    pre  = D(1,:);
    post = D(2,:);

    % Means
    Y(j,1) = mean(pre,'omitnan');
    Y(j,2) = mean(post,'omitnan');

    % Number of valid observations
    nPre  = sum(~isnan(pre));
    nPost = sum(~isnan(post));

    % SEM
    SEM(j,1) = std(pre,0,'omitnan') / sqrt(nPre);
    SEM(j,2) = std(post,0,'omitnan') / sqrt(nPost);

end

% Grouped bars
b = bar(ax,xCenters,Y,'grouped','BarWidth',0.75);

b(1).FaceColor = Cpre;
b(1).EdgeColor = 'k';
b(1).LineWidth = 0.8;

b(2).FaceColor = Cpost;
b(2).EdgeColor = 'k';
b(2).LineWidth = 0.8;

% Error bars
for k = 1:2

    x = b(k).XEndPoints;

    errorbar(ax,x,Y(:,k),SEM(:,k), ...
        'k', ...
        'LineStyle','none', ...
        'LineWidth',1.1, ...
        'CapSize',5);

end

% Separator between Inside / Outside
xline(ax,4,'--', ...
    'Color',[0.65 0.65 0.65], ...
    'LineWidth',1.2);

% Axes
xlim(ax,[0.3 7.7]);

xticks(ax,xCenters);
xticklabels(ax,objectNames);
xtickangle(ax,30);

ylabel(ax,yLabels{p}, ...
    'FontSize',13);

title(ax,panelLabels{p}, ...
    'FontWeight','bold', ...
    'FontSize',13);

set(ax, ...
    'FontSize',11, ...
    'LineWidth',1.0, ...
    'TickDir','out');

box(ax,'off');

% Add Inside / Outside labels
yl = ylim(ax);
yr = diff(yl);

ylim(ax,[yl(1)-0.14*yr yl(2)]);

yl = ylim(ax);
yText = yl(1) + 0.025*diff(yl);

text(ax,2,yText,'Inside tank', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',10);

text(ax,6,yText,'Outside tank', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',10);

% Legend only in first panel
if p == 1

    legend(ax, ...
        {'Pre-Switch','Post-Switch'}, ...
        'Location','best', ...
        'Box','off', ...
        'FontSize',10);

end
end
end % legacy combined interaction-metrics panel


%% Individual interaction figures
% plot_prepost_object_metric(freq_matrix, ...
%     'Interaction frequency (s^{-1})', 'Interaction frequency', ...
%     visibleState, figureDir, 'Interaction_Frequency_PrePost');
% 
% plot_prepost_object_metric(mean_dur_matrix, ...
%     'Mean interaction duration (s)', 'Mean interaction duration', ...
%     visibleState, figureDir, 'Interaction_MeanDuration_PrePost');
% 
% plot_prepost_object_metric(total_time_matrix, ...
%     'Total interaction time (s)', 'Total interaction time', ...
%     visibleState, figureDir, 'Interaction_TotalTime_PrePost');

frequencyChanceFile = '';
durationChanceFile = '';
totalTimeChanceFile = '';

if showSignificance
    frequencyChanceFile = fullfile( ...
        statsFolder, ...
        'InteractionFrequency_ChancePermutation.csv');

    durationChanceFile = fullfile( ...
        statsFolder, ...
        'InteractionDuration_ChancePermutation.csv');

    totalTimeChanceFile = fullfile( ...
        statsFolder, ...
        'InteractionTotalTime_ChancePermutation.csv');
end

plot_prepost_object_metric( ...
    freq_matrix, ...
    'Interaction frequency (s^{-1})', ...
    'Interaction frequency', ...
    visibleState, ...
    figureDir, ...
    'Interaction_Frequency_PrePost', ...
    frequencyChanceFile);

plot_prepost_object_metric( ...
    mean_dur_matrix, ...
    'Mean interaction duration (s)', ...
    'Mean interaction duration', ...
    visibleState, ...
    figureDir, ...
    'Interaction_MeanDuration_PrePost', ...
    durationChanceFile);

plot_prepost_object_metric( ...
    total_time_matrix, ...
    'Total interaction time (s)', ...
    'Total interaction time', ...
    visibleState, ...
    figureDir, ...
    'Interaction_TotalTime_PrePost', ...
    totalTimeChanceFile);


%% Transfer-entropy figures
plot_te_analysis(Results.TE,visibleState,figureDir);

te_metric_cells = build_te_metric_cells(Results.TE,conditionDefs);

teChanceStatsFile = '';
if showSignificance
    teChanceStatsFile = fullfile(statsFolder,'TE_ChancePermutation.csv');
end

plot_prepost_object_metric( ...
    te_metric_cells, ...
    'Transfer entropy (bits)', ...
    'Robot-to-fish transfer entropy', ...
    visibleState, ...
    figureDir, ...
    'TransferEntropy_PrePost', ...
    teChanceStatsFile);


%% Save complete MATLAB results
if ~useCachedResults
    save(cacheFile, 'Results', '-v7.3');
    fprintf('Saved results cache:\n%s\n', cacheFile);
end

fprintf('\nAnalysis/figure export complete.\n');
fprintf('Output folder: %s\n',outputFolder);
