%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [TimeSeries_Fish, TimeSeries_Robot] = Load_TimeSeries_A_and_C(Parameters,Dat)

x_max = Parameters.x_max;
y_tank = Parameters.y_tank;
y_max = Parameters.y_max;
fps = Parameters.fps;
dt = Parameters.dt;
Xmin = Parameters.Xmin;
Xmax = Parameters.Xmax;
Ymin = Parameters.Ymin;
Ymax= Parameters.Ymax;
T_total = Parameters.T_total;
n_trials = Parameters.n_trials;
scale_x = Parameters.scale_x;
scale_y = Parameters.scale_y;
x_rect = Parameters.x_rect;
y_rect = Parameters.y_rect;

px2cm_x = @(px) (2*(px - Xmin)./(Xmax-Xmin) - 1) * scale_x;
px2cm_y = @(py) (2*(py - Ymin)./(Ymax-Ymin) - 1) * scale_y;

% Fish
XX_f = nan(T_total,n_trials);
YY_f = nan(T_total,n_trials);
VV_f = nan(T_total,n_trials);
WW_f = nan(T_total,n_trials);
AA_f = nan(T_total,n_trials);
VVx_f = nan(T_total,n_trials);
VVy_f = nan(T_total,n_trials);
HH_f = nan(T_total,n_trials);

% Robot
XX_r = nan(T_total,n_trials);
YY_r = nan(T_total,n_trials);
VV_r = nan(T_total,n_trials);
WW_r = nan(T_total,n_trials);
AA_r = nan(T_total,n_trials);
VVx_r = nan(T_total,n_trials);
VVy_r = nan(T_total,n_trials);
HH_r = nan(T_total,n_trials);

NN_raw = nan(T_total,n_trials);
NN_filte = nan(T_total,n_trials);
NN_smoot = nan(T_total,n_trials);
NN_firing = nan(T_total,n_trials);

for i = 1:n_trials
    if i>length(Dat), continue; end
    Xf_px = Dat(i).Xf;
    Yf_px = Dat(i).Yf;
    Xr_px = Dat(i).Xr;
    Yr_px = Dat(i).Yr;

    n = min([T_total,length(Xf_px), length(Yf_px),length(Xr_px), length(Yr_px)]);

    Xf = px2cm_x(Xf_px(1:n));
    Yf = px2cm_y(Yf_px(1:n));
    Xr = px2cm_x(Xr_px(1:n));
    Yr = px2cm_y(Yr_px(1:n));
    N_unfil = Dat(i).Unfiltered(1:n);
    N_filte = Dat(i).Filtered(1:n);
    N_smoot = Dat(i).Smoothed(1:n);
    N_firing = Dat(i).BinarySig(1:n);

    % Filtering
    order = 4; % Filter order
    fc = 10; % Desired cutoff frequency in Hz
    Fs = 1/dt;
    wc = fc / (Fs/2); % Normalize to Nyquist
    
    [Xfilter,Yfilter,Vf,Wf,Af,Headingf,Vxf,Vyf] = Kinematic_Variables(Xf,Yf,dt,x_max/2,y_max/2,order,wc);
    [Xr_filter,Yr_filter,Vr,Wr,Ar,Headingr,Vxr,Vyr] = Kinematic_Variables(Xr,Yr,dt,x_max/2,y_max/2,order,wc);

    % save time series Fish
    XX_f(1:length(Yfilter),i) =  Xfilter;
    YY_f(1:length(Yfilter),i) =  Yfilter;
    VV_f(1:length(Vf),i) =  Vf;
    WW_f(1:length(Wf),i) =  Wf;
    AA_f(1:length(Af),i) =  Af;
    VVx_f(1:length(Vxf),i) =  Vxf;
    VVy_f(1:length(Vyf),i) =  Vyf;
    HH_f(1:length(Yfilter),i) =  Headingf;
   
    TimeSeries_Fish.XX_f = XX_f;
    TimeSeries_Fish.YY_f = YY_f;
    TimeSeries_Fish.VV_f = VV_f;
    TimeSeries_Fish.WW_f = WW_f;
    TimeSeries_Fish.AA_f = AA_f;
    TimeSeries_Fish.VVx_f = VVx_f;
    TimeSeries_Fish.VVy_f = VVy_f;
    TimeSeries_Fish.HH_f = HH_f;

    % save time series Robot
    XX_r(1:length(Xr_filter),i) =  Xr_filter;
    YY_r(1:length(Yr_filter),i) =  Yr_filter;
    VV_r(1:length(Vr),i) =  Vr;
    WW_r(1:length(Wr),i) =  Wr;
    AA_r(1:length(Ar),i) =  Ar;
    VVx_r(1:length(Vxr),i) =  Vxr;
    VVy_r(1:length(Vyr),i) =  Vyr;
    HH_r(1:length(Vyr),i) =  Headingr;

    NN_raw(1:length(N_unfil),i) = N_unfil;
    NN_filte(1:length(N_filte),i) = N_filte;
    NN_smoot(1:length(N_smoot),i) = N_smoot;
    NN_firing(1:length(N_firing),i) = N_firing;


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
end

end