% test_drift_scenarios.m
% Runs three drift scenarios and compares KPI outputs.
% Place in project root and run after verifying main.m works.

clear all; close all; clc;

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', ...
    'src/diagnostics', 'src/gates', 'src/CSACs', 'config')

config = jsondecode(fileread("config.json"));

% Shared setup
R_base = config.project.initialization.ukf.measurement_noise^2;
Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, ...
               (config.project.initialization.ukf.process_noise_U)^2]);
P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, ...
               (config.project.initialization.ukf.state_uncertainty_U)^2]);
innovation_gate_initial = config.project.initialization.ukf.state_uncertainty_offset ...
                        * config.project.initialization.innovation_gate_N_sigma;

kpi_config.offset_tolerance = 0.3;
kpi_config.U_tolerance = 0.02;
kpi_config.convergence_P_offset = 0.5;
kpi_config.convergence_P_U = 0.05;
kpi_config.convergence_hold_days = 14;

[T_soil_C, T_air_C] = soilTemp(config);
daily_T_air_max_table = build_daily_T_air_max_table(T_air_C);

%% Define scenarios
scenarios = struct();

scenarios(1).name = 'no_drift';
scenarios(1).drift_config.type = 'none';
scenarios(1).drift_config.offset_drift_per_year = 0;
scenarios(1).drift_config.step_time = NaT;
scenarios(1).drift_config.offset_step = 0;

scenarios(2).name = 'linear_drift';
scenarios(2).drift_config.type = 'linear';
scenarios(2).drift_config.offset_drift_per_year = 0.3; % 0.3 °C/year
scenarios(2).drift_config.step_time = NaT;
scenarios(2).drift_config.offset_step = 0;

scenarios(3).name = 'step_change';
scenarios(3).drift_config.type = 'step';
scenarios(3).drift_config.offset_drift_per_year = 0;
scenarios(3).drift_config.step_time = datetime(2019, 6, 1); % mid-dataset step
scenarios(3).drift_config.offset_step = 0.5; % +0.5 °C step

%% Run each scenario
networks = config.project.datasets.datasets;
% Use only the first network and first CSAC for quick testing
network_id = networks(1);
[meter_data, ~, topology] = importData(config, network_id);
timestamps = unique(meter_data.timestamp);
csac_ids = [topology.cul_de_sacs.id];
test_csac = csac_ids(1);

results = struct();

for s = 1:numel(scenarios)
    fprintf('\n========================================\n');
    fprintf('SCENARIO: %s\n', scenarios(s).name);
    fprintf('========================================\n');

    houses_csac = ([topology.houses.cul_de_sac_id] == test_csac);
    topology_csac = topology.houses(houses_csac);

    cs = initialize_csac_state(topology, topology_csac, houses_csac, ...
        meter_data, length(timestamps), Q_base, R_base, P_base, ...
        innovation_gate_initial, config, network_id);

    % Apply drift
    cs.meter_data = apply_offset_drift_to_data(...
        cs.meter_data, cs.house_ids, scenarios(s).drift_config, timestamps);

    % Generate truth
    true_traj = generate_true_trajectories(...
        cs.ground_truth, timestamps, scenarios(s).drift_config);

    % Run filter
    for t = 1:length(timestamps)
        time = timestamps(t);
        current_data = cs.meter_data(cs.meter_data.timestamp == time, :);
        current_data = sortrows(current_data, 'house_id');
        current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;

        if isempty(current_data) || isempty(current_T_soil_C)
            continue
        end

        cs = process_csac_timestep(cs, t, time, current_data, current_T_soil_C, ...
            daily_T_air_max_table, P_base, config, test_csac);
    end

    % Compute KPIs
    output_folder = fullfile('results', 'drift_test', scenarios(s).name);
    kpi_summary = compute_and_save_network_kpis(cs, test_csac, output_folder, ...
        kpi_config, true_traj);

    results(s).name = scenarios(s).name;
    results(s).kpi_summary = kpi_summary;
    results(s).cs = cs;

    fprintf('\n');
end

%% Print comparison
fprintf('\n========================================\n');
fprintf('SCENARIO COMPARISON (CSAC %d)\n', test_csac);
fprintf('========================================\n');
fprintf('%-15s | %10s | %10s | %10s | %10s\n', ...
    'Scenario', 'TW-MAE Off', 'TW-MAE U', 'Final|Off|', 'Final|U|');
fprintf('%s\n', repmat('-', 1, 65));
for s = 1:numel(results)
    k = results(s).kpi_summary;
    fprintf('%-15s | %10.4f | %10.5f | %10.4f | %10.5f\n', ...
        results(s).name, ...
        mean(k.tw_mae_offset, 'omitnan'), ...
        mean(k.tw_mae_U, 'omitnan'), ...
        mean(abs(k.final_err_offset), 'omitnan'), ...
        mean(abs(k.final_err_U), 'omitnan'));
end
fprintf('========================================\n');

fprintf('\nDrift test complete. Check results/drift_test/ for plots.\n');