%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [XX_control, YY_control, VV_control, WW_control, AA_control, VVx_control, VVy_control] = Load_TimeSeries_Control(Parameters,Dat_Control)

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


%% Control Trials
YY_control = nan(T_total,n_trials);
VV_control = nan(T_total,n_trials);
WW_control = nan(T_total,n_trials);
AA_control = nan(T_total,n_trials);
VVx_control = nan(T_total,n_trials);
VVy_control = nan(T_total,n_trials);

for i = 1:n_trials
    if i>length(Dat_Control), continue; end
    X_px = Dat_Control(i).X;
    Y_px = Dat_Control(i).Y;

    minX(i) = min(X_px);
    maxX(i) = max(X_px);
    minY(i) = min(Y_px);
    maxY(i) = max(Y_px);

    n = min([T_total,length(X_px), length(Y_px)]);
    X_cm = px2cm_x(X_px(1:n));
    Y_cm = px2cm_y(Y_px(1:n));

    order = 4; % Filter order
    fc = 10; % Desired cutoff frequency in Hz
    Fs = 1/dt;
    wc = fc / (Fs/2); % Normalize to Nyquist
    [Xfilter, Yfilter, V, W, A, heading, Vx, Vy] = Kinematic_Variables(X_cm,Y_cm,dt,x_max/2,y_max/2,order,wc);
      
    % save time series
    XX_control(1:length(Yfilter),i) =  Xfilter;
    YY_control(1:length(Yfilter),i) =  Yfilter;
    VV_control(1:length(V),i) =  V;
    WW_control(1:length(W),i) =  W;
    AA_control(1:length(A),i) =  A;
    VVx_control(1:length(Vx),i) =  Vx;
    VVy_control(1:length(Vy),i) =  Vy;
end

end