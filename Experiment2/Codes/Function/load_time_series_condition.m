%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [TimeSeries_Fish, TimeSeries_Robot] = load_time_series_condition(Parameters, Dat)

x_max = Parameters.x_max;
y_max = Parameters.y_max;
dt = Parameters.dt;
Xmin = Parameters.Xmin;
Xmax = Parameters.Xmax;
Ymin = Parameters.Ymin;
Ymax = Parameters.Ymax;
T_total = Parameters.T_total;
scale_x = Parameters.scale_x;
scale_y = Parameters.scale_y;

px2cm_x = @(px) (2*(px - Xmin)./(Xmax-Xmin) - 1) * scale_x;
px2cm_y = @(py) (2*(py - Ymin)./(Ymax-Ymin) - 1) * scale_y;

n_trials = numel(Dat);

XX_f = nan(T_total, n_trials);
YY_f = nan(T_total, n_trials);
VV_f = nan(T_total, n_trials);
WW_f = nan(T_total, n_trials);
AA_f = nan(T_total, n_trials);
VVx_f = nan(T_total, n_trials);
VVy_f = nan(T_total, n_trials);
HH_f = nan(T_total, n_trials);

XX_r = nan(T_total, n_trials);
YY_r = nan(T_total, n_trials);
VV_r = nan(T_total, n_trials);
WW_r = nan(T_total, n_trials);
AA_r = nan(T_total, n_trials);
VVx_r = nan(T_total, n_trials);
VVy_r = nan(T_total, n_trials);
HH_r = nan(T_total, n_trials);

NN_raw = nan(T_total, n_trials);
NN_filte = nan(T_total, n_trials);
NN_smoot = nan(T_total, n_trials);
NN_firing = nan(T_total, n_trials);
CmdY = nan(T_total, n_trials);
NumFrames = nan(1, n_trials);
SwitchIdx = nan(1, n_trials);
SwitchTimeSec = nan(1, n_trials);

order = 4;
fc = 10;
Fs = 1 / dt;
wc = min(fc / (Fs / 2), 0.99);

for i = 1:n_trials
    Xf_px = Dat(i).Xf(:);
    Yf_px = Dat(i).Yf(:);
    Xr_px = Dat(i).Xr(:);
    Yr_px = Dat(i).Yr(:);

    n = min([T_total, numel(Xf_px), numel(Yf_px), numel(Xr_px), numel(Yr_px)]);

    Xf = px2cm_x(Xf_px(1:n));
    Yf = px2cm_y(Yf_px(1:n));
    Xr = px2cm_x(Xr_px(1:n));
    Yr = px2cm_y(Yr_px(1:n));

    N_unfil = get_signal(Dat(i), 'Unfiltered', n);
    N_filte = get_signal(Dat(i), 'Filtered', n);
    N_smoot = get_signal(Dat(i), 'Smoothed', n);
    N_binary = get_signal(Dat(i), 'BinarySig', n);
    RobotCommandY = get_signal(Dat(i), 'RobotCommandY', n);

    [Xfilter, Yfilter, Vf, Wf, Af, Headingf, Vxf, Vyf] = ...
        kinematic_variables_local(Xf, Yf, dt, x_max/2, y_max/2, order, wc);
    [Xr_filter, Yr_filter, Vr, Wr, Ar, Headingr, Vxr, Vyr] = ...
        kinematic_variables_local(Xr, Yr, dt, x_max/2, y_max/2, order, wc);

    write_len = min(T_total, numel(Yfilter));
    XX_f(1:write_len, i) = Xfilter(1:write_len);
    YY_f(1:write_len, i) = Yfilter(1:write_len);
    HH_f(1:write_len, i) = Headingf(1:write_len);

    write_len = min(T_total, numel(Vf));
    VV_f(1:write_len, i) = Vf(1:write_len);
    AA_f(1:write_len, i) = Af(1:write_len);
    VVx_f(1:write_len, i) = Vxf(1:write_len);
    VVy_f(1:write_len, i) = Vyf(1:write_len);

    write_len = min(T_total, numel(Wf));
    WW_f(1:write_len, i) = Wf(1:write_len);

    write_len = min(T_total, numel(Xr_filter));
    XX_r(1:write_len, i) = Xr_filter(1:write_len);
    YY_r(1:write_len, i) = Yr_filter(1:write_len);
    HH_r(1:write_len, i) = Headingr(1:write_len);

    write_len = min(T_total, numel(Vr));
    VV_r(1:write_len, i) = Vr(1:write_len);
    AA_r(1:write_len, i) = Ar(1:write_len);
    VVx_r(1:write_len, i) = Vxr(1:write_len);
    VVy_r(1:write_len, i) = Vyr(1:write_len);

    write_len = min(T_total, numel(Wr));
    WW_r(1:write_len, i) = Wr(1:write_len);

    NN_raw(1:n, i) = N_unfil;
    NN_filte(1:n, i) = N_filte;
    NN_smoot(1:n, i) = N_smoot;
    NN_firing(1:n, i) = N_binary;
    CmdY(1:n, i) = RobotCommandY;
    NumFrames(i) = n;

    idxSwitch = find(RobotCommandY > -80, 1, 'first');
    if ~isempty(idxSwitch)
        SwitchIdx(i) = idxSwitch;
        SwitchTimeSec(i) = (idxSwitch - 1) * dt;
    end
end

TimeSeries_Fish.XX_f = XX_f;
TimeSeries_Fish.YY_f = YY_f;
TimeSeries_Fish.VV_f = VV_f;
TimeSeries_Fish.WW_f = WW_f;
TimeSeries_Fish.AA_f = AA_f;
TimeSeries_Fish.VVx_f = VVx_f;
TimeSeries_Fish.VVy_f = VVy_f;
TimeSeries_Fish.HH_f = HH_f;
TimeSeries_Fish.NumFrames = NumFrames;

TimeSeries_Robot.XX_r = XX_r;
TimeSeries_Robot.YY_r = YY_r;
TimeSeries_Robot.VV_r = VV_r;
TimeSeries_Robot.WW_r = WW_r;
TimeSeries_Robot.AA_r = AA_r;
TimeSeries_Robot.VVx_r = VVx_r;
TimeSeries_Robot.VVy_r = VVy_r;
TimeSeries_Robot.HH_r = HH_r;
TimeSeries_Robot.NN_raw = NN_raw;
TimeSeries_Robot.NN_filte = NN_filte;
TimeSeries_Robot.NN_smoot = NN_smoot;
TimeSeries_Robot.NN_firing = NN_firing;
TimeSeries_Robot.CmdY = CmdY;
TimeSeries_Robot.SwitchIdx = SwitchIdx;
TimeSeries_Robot.SwitchTimeSec = SwitchTimeSec;
TimeSeries_Robot.NumFrames = NumFrames;

end
