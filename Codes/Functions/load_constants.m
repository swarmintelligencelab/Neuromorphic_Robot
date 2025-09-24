%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Neuromorphic Robot Modulates Emotional Behavior in Live Fish
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Parameters] = load_constants(Dat_Control)

% Tank Dimensions;
x_max = 76;     % cm
y_tank = 50.4;  % cm 
y_max = 41;     % from bottom to water height; 

% Define Constants
fps = 55.7;
dt  = 1/fps;

T_total = 10*60*(1/dt); % Total Experimental Time, 10 min
T_five = 5*60*(1/dt);    % 5 min of experimental time

n_trials = 10; % ten trials per experiment

% calibrate tank dimnesions
for i = 1:10
    X_px = Dat_Control(i).X;
    Y_px = Dat_Control(i).Y;

    minX(i) = min(X_px);
    maxX(i) = max(X_px);
    minY(i) = min(Y_px);
    maxY(i) = max(Y_px);

end

Xmin = median(minX);
Xmax = median(maxX);
Ymin = median(minY);
Ymax = median(maxY);

scale_x = x_max/2.1;
scale_y = y_max/2.1;

% Define rectangle vertices to plot trajectories. This are the boundaries
x_rect = [-x_max/2, x_max/2, x_max/2, -x_max/2, -x_max/2]; % X-coordinates
y_rect = [-y_max/2, -y_max/2, y_max/2, y_max/2, -y_max/2]; % Y-coordinates

Parameters.x_max = x_max;
Parameters.y_tank = y_tank;
Parameters.y_max = y_max;
Parameters.fps = fps;
Parameters.dt = dt;
Parameters.Xmin = Xmin;
Parameters.Xmax = Xmax;
Parameters.Ymin = Ymin;
Parameters.Ymax = Ymax;
Parameters.scale_x = scale_x;
Parameters.scale_y= scale_y;
Parameters.x_rect = x_rect;
Parameters.y_rect = y_rect;
Parameters.T_total = T_total;
Parameters.T_five = T_five;
Parameters.n_trials = n_trials;

end