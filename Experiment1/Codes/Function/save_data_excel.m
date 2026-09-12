%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/06/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function save_data_excel(M_WW_control,M_WW_f_C2A,M_WW_f_A2C,nameFilei)
num_subjects = size(M_WW_control, 2);   % 10 subjects
num_timepoints = size(M_WW_control, 1); % 10 timepoints
time_labels = (1:num_timepoints)';      % Numeric time labels (t1, t2, ..., t10)

reshape_data = @(M, cond_name) table(...
    repelem((1:num_subjects), num_timepoints)', ... % Subject column
    repmat(time_labels, num_subjects, 1), ...       % Time column
    M(:), ...                                       % Value column
    repmat({cond_name}, num_subjects*num_timepoints, 1), ... % Condition column
    'VariableNames', {'Subject', 'Time', 'Value', 'Condition'});

T_control = reshape_data(M_WW_control, "Control");
T_C2A = reshape_data(M_WW_f_C2A, "C2A");
T_A2C = reshape_data(M_WW_f_A2C, "A2C");

long_format_table = [T_control; T_C2A; T_A2C];

writetable(long_format_table, nameFilei);

end