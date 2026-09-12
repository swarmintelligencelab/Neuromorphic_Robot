%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function RasterInfo = plot_neuron_raster_algorithms(TimeSeries_Robot, conditionDefs, Parameters, ...
    visibleState, figureDir)
windowFrames = 2 * Parameters.phaseWindowFrames;
plotTimeMin = ((0:windowFrames-1) * Parameters.dt) / 60;
nConditions = numel(TimeSeries_Robot);

rawByCondition = cell(1, nConditions);
binaryByCondition = cell(1, nConditions);

for c = 1:nConditions
    rawByCondition{c} = extract_switch_window_raster(TimeSeries_Robot{c}.NN_raw, windowFrames);
    binaryByCondition{c} = extract_switch_window_raster(TimeSeries_Robot{c}.NN_firing, windowFrames);
end

% Apply voltage threshold for raw-voltage raster
voltageThreshold = 1; % unit V

thresholdedRawByCondition = rawByCondition;

for c = 1:nConditions
    M = thresholdedRawByCondition{c};

    M(~isnan(M) & M < voltageThreshold) = 0;

    thresholdedRawByCondition{c} = M;
end

rawCLim = [0 5];
deltaByCondition = compute_raster_delta_from_pre(rawByCondition, plotTimeMin);
deltaCLim = [-0.75 0.75];
ratioByCondition = compute_raster_ratio_to_pre(rawByCondition, plotTimeMin);
ratioCLim = compute_raster_robust_clim(ratioByCondition, 2, 98);

plot_neuron_raster_pair(thresholdedRawByCondition, binaryByCondition, conditionDefs, plotTimeMin, ...
    rawCLim, 'Raw voltage', 'V_o (V)', ...
    'Switch-aligned neuron activity across six conditions', ...
    visibleState, figureDir, 'Neuron_Raster_SwitchAligned_Raw');

plot_neuron_raster_pair(deltaByCondition, binaryByCondition, conditionDefs, plotTimeMin, ...
    deltaCLim, 'Voltage change from pre-switch median', '\Delta V_o (V)', ...
    'Switch-aligned neuron activity, baseline-subtracted', ...
    visibleState, figureDir, 'Neuron_Raster_SwitchAligned_DeltaFromPre');

plot_neuron_raster_pair(ratioByCondition, binaryByCondition, conditionDefs, plotTimeMin, ...
    ratioCLim, 'Voltage ratio to pre-switch mean', 'V_o / pre mean', ...
    'Switch-aligned neuron activity, ratio to pre-switch mean', ...
    visibleState, figureDir, 'Neuron_Raster_SwitchAligned_RatioToPre');

RasterInfo = struct();
RasterInfo.rawCLim = rawCLim;
RasterInfo.deltaCLim = deltaCLim;
RasterInfo.ratioCLim = ratioCLim;
RasterInfo.windowFrames = windowFrames;
RasterInfo.windowSec = windowFrames * Parameters.dt;
RasterInfo.switchTimeSec = Parameters.phaseWindowSec;

end
