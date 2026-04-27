% test_drift_scenarios.m
% Runs drift scenarios on a single house and compares KPI outputs.
% Place in project root and run.

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

%% Define which house drifts (index within the CSAC, 1-based)
drift_house_index = 5;

%% Define scenarios — realistic drift rates and step sizes
scenarios = struct();

scenarios(1).name = 'no_drift';
scenarios(1).drift_config.type = 'none';
scenarios(1).drift_config.house_index = drift_house_index;
scenarios(1).drift_config.offset_drift_per_year = 0;
scenarios(1).drift_config.step_time = NaT;
scenarios(1).drift_config.offset_step = 0;

% Slow drift — barely detectable, realistic ageing
scenarios(2).name = 'linear_0.05C_yr';
scenarios(2).drift_config.type = 'linear';
scenarios(2).drift_config.house_index = drift_house_index;
scenarios(2).drift_config.offset_drift_per_year = 0.05;
scenarios(2).drift_config.step_time = NaT;
scenarios(2).drift_config.offset_step = 0;

% Moderate drift — detectable over a few seasons
scenarios(3).name = 'linear_0.1C_yr';
scenarios(3).drift_config.type = 'linear';
scenarios(3).drift_config.house_index = drift_house_index;
scenarios(3).drift_config.offset_drift_per_year = 0.1;
scenarios(3).drift_config.step_time = NaT;
scenarios(3).drift_config.offset_step = 0;

% Faster drift — should be clearly visible
scenarios(4).name = 'linear_0.3C_yr';
scenarios(4).drift_config.type = 'linear';
scenarios(4).drift_config.house_index = drift_house_index;
scenarios(4).drift_config.offset_drift_per_year = 0.3;
scenarios(4).drift_config.step_time = NaT;
scenarios(4).drift_config.offset_step = 0;

% Small step — partial sensor degradation
scenarios(5).name = 'step_0.3C';
scenarios(5).drift_config.type = 'step';
scenarios(5).drift_config.house_index = drift_house_index;
scenarios(5).drift_config.offset_drift_per_year = 0;
scenarios(5).drift_config.step_time = datetime(2019, 1, 15);
scenarios(5).drift_config.offset_step = 0.3;

% Medium step — noticeable sensor issue
scenarios(6).name = 'step_0.5C';
scenarios(6).drift_config.type = 'step';
scenarios(6).drift_config.house_index = drift_house_index;
scenarios(6).drift_config.offset_drift_per_year = 0;
scenarios(6).drift_config.step_time = datetime(2019, 1, 15);
scenarios(6).drift_config.offset_step = 0.5;

% Large step — sensor replacement or failure
scenarios(7).name = 'step_1.0C';
scenarios(7).drift_config.type = 'step';
scenarios(7).drift_config.house_index = drift_house_index;
scenarios(7).drift_config.offset_drift_per_year = 0;
scenarios(7).drift_config.step_time = datetime(2019, 1, 15);
scenarios(7).drift_config.offset_step = 1.0;

%% Run each scenario
networks = config.project.datasets.datasets;
network_id = networks(1);
[meter_data, ~, topology] = importData(config, network_id);
timestamps = unique(meter_data.timestamp);
csac_ids = [topology.cul_de_sacs.id];
test_csac = csac_ids(1);

results = struct();

for s = 1:numel(scenarios)
    fprintf('\n========================================\n');
    fprintf('SCENARIO %d/%d: %s\n', s, numel(scenarios), scenarios(s).name);
    fprintf('========================================\n');

    houses_csac = ([topology.houses.cul_de_sac_id] == test_csac);
    topology_csac = topology.houses(houses_csac);

    cs = initialize_csac_state(topology, topology_csac, houses_csac, ...
        meter_data, length(timestamps), Q_base, R_base, P_base, ...
        innovation_gate_initial, config, network_id);

    cs.meter_data = apply_offset_drift_to_data(...
        cs.meter_data, cs.house_ids, scenarios(s).drift_config, timestamps);

    true_traj = generate_true_trajectories(...
        cs.ground_truth, timestamps, scenarios(s).drift_config);

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

    output_folder = fullfile('results', 'drift_test', scenarios(s).name);
    kpi_summary = compute_and_save_network_kpis(cs, test_csac, output_folder, ...
        kpi_config, true_traj);

    results(s).name = scenarios(s).name;
    results(s).kpi_summary = kpi_summary;
    results(s).cs = cs;
    results(s).true_traj = true_traj;
end

%% Print comparison — drifting house only
di = drift_house_index;
fprintf('\n========================================\n');
fprintf('DRIFT HOUSE COMPARISON (house index %d, CSAC %d)\n', di, test_csac);
fprintf('========================================\n');
fprintf('%-18s | %10s | %10s | %10s | %10s | %10s\n', ...
    'Scenario', 'TW-MAE', 'Final|Err|', 'SS-MAE', 'MaxExc', 'Stability');
fprintf('%s\n', repmat('-', 1, 78));
for s = 1:numel(results)
    k = results(s).kpi_summary;
    fprintf('%-18s | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        results(s).name, ...
        k.tw_mae_offset(di), ...
        abs(k.final_err_offset(di)), ...
        k.ss_mae_offset(di), ...
        k.max_exc_offset(di), ...
        k.stability_offset(di));
end
fprintf('========================================\n');

%% Print comparison — non-drifting houses (should be unaffected)
other_idx = setdiff(1:size(results(1).kpi_summary, 1), di);
fprintf('\nNON-DRIFTING HOUSES (mean across %d houses)\n', numel(other_idx));
fprintf('%-18s | %10s | %10s | %10s\n', 'Scenario', 'TW-MAE', 'Final|Err|', 'Stability');
fprintf('%s\n', repmat('-', 1, 55));
for s = 1:numel(results)
    k = results(s).kpi_summary;
    fprintf('%-18s | %10.4f | %10.4f | %10.4f\n', ...
        results(s).name, ...
        mean(k.tw_mae_offset(other_idx), 'omitnan'), ...
        mean(abs(k.final_err_offset(other_idx)), 'omitnan'), ...
        mean(k.stability_offset(other_idx), 'omitnan'));
end
fprintf('========================================\n');

%% Plot trajectory of drifting house across scenarios
fig = figure('Visible', 'on', 'Position', [100, 100, 1400, 700]);
colors = lines(numel(scenarios));

for s = 1:numel(scenarios)
    cs = results(s).cs;
    est_offset = squeeze(cs.logger.state_estimates(1, di, :));
    true_offset = results(s).true_traj{di}.offset;
    est_std = squeeze(sqrt(cs.logger.covariance_posterior(1, di, :)));

    valid = ~isnat(cs.logger.timestamps);
    t_plot = cs.logger.timestamps(valid);
    est_plot = est_offset(valid);
    true_plot = true_offset(valid);
    std_plot = est_std(valid);

    subplot(3, 1, 1);
    hold on;
    plot(t_plot, est_plot, '-', 'Color', colors(s,:), 'LineWidth', 1.2, ...
        'DisplayName', results(s).name);

    subplot(3, 1, 2);
    hold on;
    plot(t_plot, true_plot, '-', 'Color', colors(s,:), 'LineWidth', 1.2, ...
        'DisplayName', results(s).name);

    subplot(3, 1, 3);
    hold on;
    err_plot = est_plot(:) - true_plot(:);
    plot(t_plot, err_plot, '-', 'Color', colors(s,:), 'LineWidth', 1.2, ...
        'DisplayName', results(s).name);
end

subplot(3, 1, 1);
ylabel('Estimated Offset [°C]');
title(sprintf('House %d (index %d) — Estimated Offset', cs.house_ids(di), di));
legend('Location', 'best', 'FontSize', 7);
grid on;

subplot(3, 1, 2);
ylabel('True Offset [°C]');
title('True Offset Trajectory');
legend('Location', 'best', 'FontSize', 7);
grid on;

subplot(3, 1, 3);
ylabel('Error [°C]');
xlabel('Time');
title('Estimation Error (Estimated - True)');
legend('Location', 'best', 'FontSize', 7);
yline(0, 'k--', 'HandleVisibility', 'off');
grid on;

output_folder = fullfile('results', 'drift_test');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
saveas(fig, fullfile(output_folder, 'drift_trajectory_comparison.png'));

fprintf('\nDrift test complete. Check results/drift_test/ for plots.\n');