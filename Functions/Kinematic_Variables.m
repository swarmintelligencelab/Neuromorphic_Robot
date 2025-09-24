%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Xfilter,Yfilter,V,W,A,heading,Vx,Vy] = Kinematic_Variables(Xi,Yi,dt,Xmax,Ymax,order,wc)

% Outlier removal
speed_thresh = 100;        % cm/s
accel_thresh = 1000;       % cm/s²

vx_raw = gradient(Xi, dt);
vy_raw = gradient(Yi, dt);
V_raw = sqrt(vx_raw.^2 + vy_raw.^2);
A_raw = gradient(V_raw, dt);  % Acceleration magnitude (scalar)

outlier_speed = V_raw > speed_thresh;
outlier_accel = abs(A_raw) > accel_thresh;
outlier_idx = outlier_speed | outlier_accel;
Xi(outlier_idx) = NaN;
Yi(outlier_idx) = NaN;

% Apply spatial bounds
Xi(abs(Xi) > Xmax) = NaN;
Yi(abs(Yi) > Ymax) = NaN;

% Count removals
n_total = length(Xi);
n_combined = sum(isnan(Xi));

X = fillmissing(Xi, 'linear');
Y = fillmissing(Yi, 'linear');

[b, a] = butter(order, wc, 'low');
Xfilter = filtfilt(b, a, X);
Yfilter = filtfilt(b, a, Y);

% Cap the values within [-Xmax, Xmax] and [-Ymax, Ymax]
Xfilter = max(min(Xfilter, Xmax), -Xmax);
Yfilter = max(min(Yfilter, Ymax), -Ymax);

Vx = gradient(Xfilter,dt);
Vy = gradient(Yfilter,dt);
V = sqrt(Vx.^2 + Vy.^2);

A = gradient(V,dt);

heading = wrapToPi(atan2(Vy, Vx)); % heading computation
W = turnrate_comp(cos(heading),sin(heading),dt);

end