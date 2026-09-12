%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Xfilter, Yfilter, V, W, A, heading, Vx, Vy] = ...
    kinematic_variables_local(Xi, Yi, dt, Xmax, Ymax, order, wc)

speed_thresh = 100;
accel_thresh = 1000;

vx_raw = gradient(Xi, dt);
vy_raw = gradient(Yi, dt);
V_raw = hypot(vx_raw, vy_raw);
A_raw = gradient(V_raw, dt);

outlier_idx = V_raw > speed_thresh | abs(A_raw) > accel_thresh;
Xi(outlier_idx) = NaN;
Yi(outlier_idx) = NaN;

Xi(abs(Xi) > Xmax) = NaN;
Yi(abs(Yi) > Ymax) = NaN;

if all(isnan(Xi)), Xi(:) = 0; end
if all(isnan(Yi)), Yi(:) = 0; end

X = fillmissing(Xi, 'linear', 'EndValues', 'nearest');
Y = fillmissing(Yi, 'linear', 'EndValues', 'nearest');

[b, a] = butter(order, wc, 'low');
Xfilter = filtfilt(b, a, X);
Yfilter = filtfilt(b, a, Y);

Xfilter = max(min(Xfilter, Xmax), -Xmax);
Yfilter = max(min(Yfilter, Ymax), -Ymax);

Vx = gradient(Xfilter, dt);
Vy = gradient(Yfilter, dt);
V = hypot(Vx, Vy);
A = gradient(V, dt);

heading = atan2(Vy, Vx);
W = turnrate_local(heading, dt);

end
