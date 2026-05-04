% test_reference_sensor.m
% Tests the impact of reference sensor placement on estimation accuracy.
% Runs over multiple networks and tests single + dual sensor placements.

clear all; close all; clc;

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', ...
    'src/diagnostics', 'src/gates', 'src/CSACs', 'config')

config = jsondecode(fileread("config.json"));

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

drift_config.type = 'none';
drift_config.house_index = 1;
drift_config.offset_drift_per_year = 0;
drift_config.step_time = NaT;
drift_config.offset_step = 0;

[T_soil_C, T_air_C] = soilTemp(config);
daily_T_air_max_table = build_daily_T_air_max_table(T_air_C);

measurement_noise = config.project.initialization.ukf.measurement_noise;
sensor_noise = 0.1;

networks_to_test = config.project.datasets.datasets(:)';
num_networks = numel(networks_to_test);

%% Build scenario list
% We need: no sensor, 4 single sensors, 6 dual sensor combinations
% CSAC IDs will be determined per network, but we define scenarios by index
scenario_defs = struct();
idx = 1;

% No sensor
scenario_defs(idx).name = 'no_sensor';
scenario_defs(idx).sensor_indices = [];
idx = idx + 1;

% Single sensors (indices 1-4 for CSACs)
for c = 1:4
    scenario_defs(idx).name = sprintf('single_csac_%d', c-1);
    scenario_defs(idx).sensor_indices = c;
    idx = idx + 1;
end

% Dual sensor combinations
combos = nchoosek(1:4, 2);
for k = 1:size(combos, 1)
    scenario_defs(idx).name = sprintf('dual_csac_%d_%d', combos(k,1)-1, combos(k,2)-1);
    scenario_defs(idx).sensor_indices = combos(k,:);
    idx = idx + 1;
end

num_scenarios = numel(scenario_defs);
fprintf('Testing %d scenarios across %d networks\n', num_scenarios, num_networks);

%% Pre-allocate aggregate results
agg_results = struct();
for s = 1:num_scenarios
    agg_results(s).name = scenario_defs(s).name;
    agg_results(s).all_tw_mae_offset = [];
    agg_results(s).all_final_err_offset = [];
    agg_results(s).all_tw_mae_U = [];
    agg_results(s).all_final_err_U = [];
    agg_results(s).all_offset_bias = [];
    agg_results(s).all_offset_std = [];
    agg_results(s).all_U_csac_err = [];
    agg_results(s).all_rej_pct = [];
    agg_results(s).all_per_csac_final_err = [];
    agg_results(s).all_per_csac_bias = [];
end

%% Main loop over networks
for net_idx = 1:num_networks
    network_id = networks_to_test(net_idx);
    fprintf('\n############################################################\n');
    fprintf('NETWORK %d (%d/%d)\n', network_id, net_idx, num_networks);
    fprintf('############################################################\n');

    [meter_data, network_data, topology] = importData(config, network_id);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];
    num_csacs = numel(csac_ids);

    U_csac_true = topology.pipe_parameters.csac_pipe.insulation_W_m_K;
    U_main_true = topology.pipe_parameters.main_pipe.insulation_W_m_K;

    % Precompute junction positions
    junction_positions = zeros(num_csacs, 1);
    for c = 1:num_csacs
        topo_idx = find([topology.cul_de_sacs.id] == csac_ids(c));
        junction_positions(c) = topology.cul_de_sacs(topo_idx).dist_on_main_m;
    end

    % Load all sensor data columns at once
    all_sensor_data = cell(num_csacs, 1);
    for c = 1:num_csacs
        col_name = sprintf('J_Main_%d_s', csac_ids(c));
        if ismember(col_name, network_data.Properties.VariableNames)
            rng(network_id * 10000 + 5e6 + csac_ids(c), 'twister');
            true_temps = network_data.(col_name);
            noise = sensor_noise * randn(size(true_temps));
            all_sensor_data{c} = struct('timestamps', network_data.timestamp, ...
                'temperatures', true_temps + noise, 'valid', true);
        else
            all_sensor_data{c} = struct('timestamps', [], 'temperatures', [], 'valid', false);
        end
    end
    rng('shuffle');

    %% Run each scenario
    for s = 1:num_scenarios
        fprintf('\n--- Network %d, Scenario %d/%d: %s ---\n', ...
            network_id, s, num_scenarios, scenario_defs(s).name);

        % Initialize all CSACs
        [all_cs, all_true_traj] = initialize_all_csacs(topology, meter_data, timestamps, ...
            Q_base, R_base, P_base, innovation_gate_initial, config, network_id, drift_config);

        shared_U_csac = all_cs{1}.U_csac;
        for c = 1:num_csacs
            all_cs{c}.U_csac = shared_U_csac;
        end

        shared_U_main = U_main_true;
        U_main_uncertainty = 0.01;

        U_csac_cfg = config.project.csac_U_estimation;
        U_csac_update_counter = 0;

        sensor_indices = scenario_defs(s).sensor_indices;
        has_sensors = ~isempty(sensor_indices);

        % Timestep loop
        for t = 1:length(timestamps)
            time = timestamps(t);
            current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
            if isempty(current_T_soil_C)
                continue
            end

            % Step 1: Process each CSAC independently
            any_csac_active = false;
            for c = 1:num_csacs
                cs = all_cs{c};
                cs.T_inlet_from_main = NaN;

                current_data = cs.meter_data(cs.meter_data.timestamp == time, :);
                current_data = sortrows(current_data, 'house_id');
                if isempty(current_data)
                    continue
                end

                cs = process_csac_timestep(cs, t, time, current_data, current_T_soil_C, ...
                    daily_T_air_max_table, P_base, config, csac_ids(c));
                all_cs{c} = cs;

                if cs.season_state.active
                    any_csac_active = true;
                end
            end

            if ~any_csac_active
                continue
            end

            % Steps 2-4: Sensor fusion (if sensors available)
            if has_sensors
                % Get flows
                junction_flows = zeros(num_csacs, 1);
                for c = 1:num_csacs
                    junction_flows(c) = all_cs{c}.current_total_flow;
                end

                % For each sensor, propagate and accumulate weighted estimates
                T_sensor_combined = nan(num_csacs, 1);
                W_sensor_combined = zeros(num_csacs, 1);

                for si = 1:numel(sensor_indices)
                    sensor_c = sensor_indices(si);
                    sensor_pos = junction_positions(sensor_c);

                    % Get sensor reading at this timestep
                    ref_T = NaN;
                    if all_sensor_data{sensor_c}.valid
                        ref_idx = find(all_sensor_data{sensor_c}.timestamps == time, 1);
                        if ~isempty(ref_idx)
                            ref_T = all_sensor_data{sensor_c}.temperatures(ref_idx);
                        end
                    end

                    if ~isfinite(ref_T)
                        continue
                    end

                    % Propagate this sensor
                    [T_prop, T_unc] = propagate_sensor_temperature(...
                        ref_T, sensor_pos, sensor_noise, ...
                        junction_positions, junction_flows, ...
                        shared_U_main, U_main_uncertainty, current_T_soil_C, topology);

                    % Accumulate inverse-variance weighted contributions
                    for c = 1:num_csacs
                        if isfinite(T_prop(c)) && isfinite(T_unc(c)) && T_unc(c) > 0
                            w = 1 / T_unc(c)^2;
                            if isnan(T_sensor_combined(c))
                                T_sensor_combined(c) = w * T_prop(c);
                                W_sensor_combined(c) = w;
                            else
                                T_sensor_combined(c) = T_sensor_combined(c) + w * T_prop(c);
                                W_sensor_combined(c) = W_sensor_combined(c) + w;
                            end
                        end
                    end
                end

                % Finalize combined sensor estimate
                for c = 1:num_csacs
                    if W_sensor_combined(c) > 0
                        T_sensor_final = T_sensor_combined(c) / W_sensor_combined(c);
                        sigma_sensor_final = 1 / sqrt(W_sensor_combined(c));

                        if isfield(all_cs{c}, 'T_inlet_fitted') && isfinite(all_cs{c}.T_inlet_fitted)
                            T_own = all_cs{c}.T_inlet_fitted;
                            sigma_csac = all_cs{c}.T_inlet_sigma;
                            sigma_csac_total = sqrt(sigma_csac^2 + ...
                                (measurement_noise / sqrt(max(1, all_cs{c}.num_houses)))^2);
                            sigma_csac_total = max(sigma_csac_total, 0.01);
                            sigma_sensor_final = max(sigma_sensor_final, 0.01);

                            w_csac = 1 / sigma_csac_total^2;
                            w_sensor = 1 / sigma_sensor_final^2;
                            T_merged = (w_csac * T_own + w_sensor * T_sensor_final) / ...
                                (w_csac + w_sensor);

                            all_cs{c}.T_inlet_from_main = T_merged;
                        end
                    end
                end

                % Step 4: Re-process CSACs with merged T_inlet
                for c = 1:num_csacs
                    cs = all_cs{c};
                    if ~isfield(cs, 'T_inlet_from_main') || ~isfinite(cs.T_inlet_from_main)
                        continue
                    end

                    current_data = cs.meter_data(cs.meter_data.timestamp == time, :);
                    current_data = sortrows(current_data, 'house_id');
                    if isempty(current_data)
                        continue
                    end

                    cs = process_csac_timestep(cs, t, time, current_data, current_T_soil_C, ...
                        daily_T_air_max_table, P_base, config, csac_ids(c));
                    all_cs{c} = cs;
                end
            end

            % U_csac estimation
            if U_csac_cfg.enabled
                U_csac_update_counter = U_csac_update_counter + 1;
                if U_csac_update_counter >= U_csac_cfg.warmup_timesteps && ...
                   mod(U_csac_update_counter, U_csac_cfg.update_interval_timesteps) == 0
                    [U_csac_new, ~] = update_shared_U_csac(all_cs, shared_U_csac, config);
                    shared_U_csac = U_csac_new;
                    for c = 1:num_csacs
                        all_cs{c}.U_csac = shared_U_csac;
                    end
                end
            end
        end

        %% Compute KPIs
        all_kpi = table();
        for c = 1:num_csacs
            kpi = compute_and_save_network_kpis(all_cs{c}, csac_ids(c), ...
                fullfile('results', 'sensor_test', sprintf('net%d', network_id), ...
                scenario_defs(s).name), kpi_config, all_true_traj{c}, network_id);
            all_kpi = [all_kpi; kpi]; %#ok<AGROW>
        end

        % Compute additional metrics
        all_offsets = [];
        all_true_offsets = [];
        for c = 1:num_csacs
            cs = all_cs{c};
            for i = 1:cs.num_houses
                all_offsets(end+1) = cs.ukf_states{i}.x(1); %#ok<AGROW>
                all_true_offsets(end+1) = cs.ground_truth.true_offset(i); %#ok<AGROW>
            end
        end
        offset_errors = all_offsets - all_true_offsets;

        total_accept = 0; total_reject = 0;
        per_csac_final_err = nan(1, num_csacs);
        per_csac_bias = nan(1, num_csacs);
        for c = 1:num_csacs
            total_accept = total_accept + all_cs{c}.gate_accept_count;
            total_reject = total_reject + all_cs{c}.gate_reject_count;
            csac_mask = all_kpi.csac_id == csac_ids(c);
            per_csac_final_err(c) = mean(abs(all_kpi.final_err_offset(csac_mask)), 'omitnan');
            cs = all_cs{c};
            offs = zeros(cs.num_houses, 1);
            for i = 1:cs.num_houses
                offs(i) = cs.ukf_states{i}.x(1);
            end
            per_csac_bias(c) = mean(offs - cs.ground_truth.true_offset);
        end
        rej_pct = 100 * total_reject / max(1, total_accept + total_reject);

        % Store aggregate results
        agg_results(s).all_tw_mae_offset(end+1) = mean(all_kpi.tw_mae_offset, 'omitnan');
        agg_results(s).all_final_err_offset(end+1) = mean(abs(all_kpi.final_err_offset), 'omitnan');
        agg_results(s).all_tw_mae_U(end+1) = mean(all_kpi.tw_mae_U, 'omitnan');
        agg_results(s).all_final_err_U(end+1) = mean(abs(all_kpi.final_err_U), 'omitnan');
        agg_results(s).all_offset_bias(end+1) = mean(offset_errors);
        agg_results(s).all_offset_std(end+1) = std(offset_errors);
        agg_results(s).all_U_csac_err(end+1) = shared_U_csac - U_csac_true;
        agg_results(s).all_rej_pct(end+1) = rej_pct;
        agg_results(s).all_per_csac_final_err = [agg_results(s).all_per_csac_final_err; per_csac_final_err];
        agg_results(s).all_per_csac_bias = [agg_results(s).all_per_csac_bias; per_csac_bias];

        fprintf('  Off MAE=%.3f, Final=%.3f, Bias=%+.3f, Uc_err=%+.4f, Rej=%.1f%%\n', ...
            mean(all_kpi.tw_mae_offset, 'omitnan'), ...
            mean(abs(all_kpi.final_err_offset), 'omitnan'), ...
            mean(offset_errors), ...
            shared_U_csac - U_csac_true, rej_pct);
    end
end

%% Print aggregate comparison table
fprintf('\n\n############################################################\n');
fprintf('AGGREGATE RESULTS ACROSS %d NETWORKS\n', num_networks);
fprintf('############################################################\n\n');

fprintf('%-20s | %8s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'Scenario', 'Off MAE', 'Off Fin', 'Off Bias', 'Off Std', 'U MAE', 'Uc err', 'Rej %%');
fprintf('%s\n', repmat('-', 1, 100));

for s = 1:num_scenarios
    a = agg_results(s);
    fprintf('%-20s | %8.3f | %8.3f | %+8.3f | %8.3f | %8.4f | %+8.4f | %7.1f\n', ...
        a.name, ...
        mean(a.all_tw_mae_offset), ...
        mean(a.all_final_err_offset), ...
        mean(a.all_offset_bias), ...
        mean(a.all_offset_std), ...
        mean(a.all_tw_mae_U), ...
        mean(a.all_U_csac_err), ...
        mean(a.all_rej_pct));
end
fprintf('========================================\n');

%% Per-CSAC breakdown (averaged across networks)
fprintf('\nPER-CSAC FINAL OFFSET ERROR (averaged across %d networks)\n', num_networks);
fprintf('%-20s | %8s | %8s | %8s | %8s\n', 'Scenario', 'CSAC0', 'CSAC1', 'CSAC2', 'CSAC3');
fprintf('%s\n', repmat('-', 1, 60));
for s = 1:num_scenarios
    a = agg_results(s);
    avg_per_csac = mean(a.all_per_csac_final_err, 1);
    fprintf('%-20s | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
        a.name, avg_per_csac(1), avg_per_csac(2), avg_per_csac(3), avg_per_csac(4));
end
fprintf('========================================\n');

%% Improvement summary
fprintf('\nIMPROVEMENT vs NO SENSOR (Final offset error)\n');
fprintf('%-20s | %10s | %10s\n', 'Scenario', 'Final Err', 'Improvement');
fprintf('%s\n', repmat('-', 1, 45));
baseline = mean(agg_results(1).all_final_err_offset);
for s = 1:num_scenarios
    fin = mean(agg_results(s).all_final_err_offset);
    improvement = 100 * (baseline - fin) / baseline;
    fprintf('%-20s | %10.3f | %+9.1f%%\n', agg_results(s).name, fin, improvement);
end
fprintf('========================================\n');

%% Save results to CSV
output_folder = fullfile('results', 'sensor_test');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

summary = table();
for s = 1:num_scenarios
    a = agg_results(s);
    row = table();
    row.scenario = {a.name};
    row.num_sensors = numel(scenario_defs(s).sensor_indices);
    row.mean_tw_mae_offset = mean(a.all_tw_mae_offset);
    row.mean_final_err_offset = mean(a.all_final_err_offset);
    row.mean_offset_bias = mean(a.all_offset_bias);
    row.mean_offset_std = mean(a.all_offset_std);
    row.mean_tw_mae_U = mean(a.all_tw_mae_U);
    row.mean_U_csac_err = mean(a.all_U_csac_err);
    row.mean_rej_pct = mean(a.all_rej_pct);
    row.num_networks = num_networks;
    summary = [summary; row]; %#ok<AGROW>
end
writetable(summary, fullfile(output_folder, 'sensor_comparison_aggregate.csv'));

fprintf('\nResults saved to %s\n', fullfile(output_folder, 'sensor_comparison_aggregate.csv'));
fprintf('Test complete.\n');