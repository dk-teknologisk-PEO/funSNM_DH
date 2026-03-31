
% inspect_data.m
%
% This script loads data for a specific CSAC, processes physical variables
% for ALL houses within it, and saves the output to a CSV file. It also
% generates a diagnostic plot for the first house.

clear all;
close all;
clc;

% Add paths to your source functions
addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'config')

% --- Configuration ---
config = jsondecode(fileread("config.json"));
t_start = config.project.time.start;
t_end = config.project.time.end;
network_to_inspect = 1;
csac_to_inspect = 1;
output_folder = fullfile('results', 'inspection_data'); % Define an output folder
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% --- Load Data ---
fprintf('Loading data for Network %d...\n', network_to_inspect);
[meter_data, network_data, topology] = importData(config, network_to_inspect);
[T_soil_C, T_air_C] = soilTemp(config);

timestamps = unique(meter_data.timestamp);

T_soil_C = T_soil_C((T_soil_C.time>=min(timestamps)) & (T_soil_C.time<=max(timestamps)),:);
T_air_C = T_air_C((T_air_C.time>=min(timestamps)) & (T_air_C.time<=max(timestamps)),:);
% Isolate the specific CSAC
houses_csac_indices = ([topology.houses.cul_de_sac_id] == csac_to_inspect);
topology_csac = topology.houses(houses_csac_indices);
houses_csac_ids = [topology_csac.id];

% --- Process Data for ALL Houses in the CSAC ---
fprintf('Processing data for all %d houses in CSAC %d...\n', length(houses_csac_ids), csac_to_inspect);
all_houses_data = []; % Initialize an empty table to append to

% Join temperature data once
T_air_C_timed = table(T_air_C.time, T_air_C.values, 'VariableNames', {'timestamp', 'T_air_C'});
T_soil_C_timed = table(T_soil_C.time, T_soil_C.values, 'VariableNames', {'timestamp', 'T_soil_C'});
timeline_data = table(timestamps, 'VariableNames', {'timestamp'});
timeline_data = outerjoin(timeline_data, T_air_C_timed, 'Keys', 'timestamp', 'MergeKeys', true);
timeline_data = outerjoin(timeline_data, T_soil_C_timed, 'Keys', 'timestamp', 'MergeKeys', true);
timeline_data.T_air_C = fillmissing(timeline_data.T_air_C, 'previous');
timeline_data.T_soil_C = fillmissing(timeline_data.T_soil_C, 'previous');

for i = 1:length(houses_csac_ids)
    house_id = houses_csac_ids(i);
    
    % Filter meter data for this specific house
    house_meter_data = meter_data(meter_data.house_id == house_id, :);
    
    % Create a complete timeline for the house
    house_timeline = outerjoin(timeline_data, house_meter_data, 'Keys', 'timestamp', 'MergeKeys', true);
    house_timeline.flow_kg_h(isnan(house_timeline.flow_kg_h)) = 0; % Assume no data means zero flow
    
    % Estimate main pipe temperature and Delta_T
    U_guess = topology_csac(i).service_pipe_insulation_W_m_K; % A typical value for a service pipe
    L_pipe = topology_csac(i).service_pipe_len_m;
    house_timeline.T_main_est_C = get_main_temp(house_timeline.T_supply_C, house_timeline.flow_kg_h, U_guess, L_pipe, house_timeline.T_soil_C);
    house_timeline.Delta_T_C = house_timeline.T_main_est_C - house_timeline.T_soil_C;

    % Append to the master table
    all_houses_data = [all_houses_data; house_timeline];
end
all_houses_data(isnan(all_houses_data.house_id),:)=[];

% --- Save the Comprehensive Data to CSV ---
output_filename = fullfile(output_folder, sprintf('physical_inputs_network%d_csac%d.csv', network_to_inspect, csac_to_inspect));
writetable(all_houses_data, output_filename);
fprintf('Saved all physical input data to %s\n', output_filename);

for i=1:length(houses_csac_ids)
    % --- Generate the Diagnostic Plot for the First House (as before) ---
    fprintf('Generating diagnostic plot for first house (House %d)...\n', houses_csac_ids(i));
    plot_data = all_houses_data(all_houses_data.house_id == houses_csac_ids(i), :);
    
    fig = figure('Name', 'Input Data Diagnostics', 'Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('Diagnostic Inputs for House %d (Network %d, CSAC %d)', houses_csac_ids(i), network_to_inspect, csac_to_inspect));
    
    % Plot 1: Flow Rate
    ax1 = subplot(3, 1, 1);
    plot(plot_data.timestamp, plot_data.flow_kg_h, 'b-');
    hold on;
    yline(config.project.cutoff.flow_cutoff, 'r--', 'Label', 'Flow Cutoff');
    title('Flow Rate');
    ylabel('Flow [kg/h]');
    grid on;
    legend('Flow Rate', 'Gating Threshold');
    
    % Plot 2: Key Temperatures
    ax2 = subplot(3, 1, 2);
    plot(plot_data.timestamp, plot_data.T_air_C, 'g-', 'DisplayName', 'Air Temperature');
    hold on;
    plot(plot_data.timestamp, plot_data.T_soil_C, 'Color', [0.5 0.2 0], 'DisplayName', 'Soil Temperature');
    yline(config.project.initialization.max_air_temperature, 'r--', 'Label', sprintf('Air Temp Cutoff (%.0f°C)', config.project.initialization.max_air_temperature));
    title('Ambient Temperatures');
    ylabel('Temperature [°C]');
    grid on;
    legend('show');
    
    % Plot 3: The Driving Delta-T
    ax3 = subplot(3, 1, 3);
    plot(plot_data.timestamp, plot_data.Delta_T_C, 'm-');
    hold on;
    yline(config.project.cutoff.delta_T_gate_threshold, 'r--', 'Label', 'Delta-T Cutoff');
    title('Driving Temperature Difference (T_{main} - T_{soil})');
    ylabel('\Delta T [°C]');
    xlabel('Date');
    grid on;
    legend('Estimated \Delta T', 'Gating Threshold');
    
    % Link axes for synchronized zooming and panning
    linkaxes([ax1, ax2, ax3], 'x');
    xlim([datetime(t_start), datetime(t_end)]);
    output_filename = fullfile(output_folder, sprintf('physical_inputs_house%d_network%d_csac%d.png', houses_csac_ids(i), network_to_inspect, csac_to_inspect));
    saveas(fig, output_filename);
    close all
end
disp('Plot generated. CSV files are ready for analysis.');
