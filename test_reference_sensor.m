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
output_folder = fullfile('results', 'sensor_test');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% Long-form table with one row per (network, scenario) for global strategy analysis.
all_scenario_rows = table();

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

    % Build scenarios dynamically for this network size.
    scenario_defs = struct();
    idx = 1;

    scenario_defs(idx).name = 'no_sensor';
    scenario_defs(idx).sensor_indices = [];
    idx = idx + 1;

    for c = 1:num_csacs
        scenario_defs(idx).name = sprintf('single_csac_%d', csac_ids(c));
        scenario_defs(idx).sensor_indices = c;
        idx = idx + 1;
    end

    if num_csacs >= 2
        combos = nchoosek(1:num_csacs, 2);
        for k = 1:size(combos, 1)
            c1 = csac_ids(combos(k, 1));
            c2 = csac_ids(combos(k, 2));
            scenario_defs(idx).name = sprintf('dual_csac_%d_%d', c1, c2);
            scenario_defs(idx).sensor_indices = combos(k, :);
            idx = idx + 1;
        end
    end

    num_scenarios = numel(scenario_defs);
    fprintf('Testing %d scenarios for network %d (num_csacs=%d)\n', ...
        num_scenarios, network_id, num_csacs);

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

    net_results = struct();
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

        % Store per-network results
        net_results(s).name = scenario_defs(s).name;
        net_results(s).num_sensors = numel(scenario_defs(s).sensor_indices);
        net_results(s).tw_mae_offset = mean(all_kpi.tw_mae_offset, 'omitnan');
        net_results(s).final_err_offset = mean(abs(all_kpi.final_err_offset), 'omitnan');
        net_results(s).offset_bias = mean(offset_errors);
        net_results(s).offset_std = std(offset_errors);
        net_results(s).tw_mae_U = mean(all_kpi.tw_mae_U, 'omitnan');
        net_results(s).final_err_U = mean(abs(all_kpi.final_err_U), 'omitnan');
        net_results(s).U_csac_err = shared_U_csac - U_csac_true;
        net_results(s).rej_pct = rej_pct;
        net_results(s).per_csac_final_err = per_csac_final_err;
        net_results(s).per_csac_bias = per_csac_bias;
        net_results(s).csac_ids = csac_ids;

        fprintf('  Off MAE=%.3f, Final=%.3f, Bias=%+.3f, Uc_err=%+.4f, Rej=%.1f%%\n', ...
            net_results(s).tw_mae_offset, net_results(s).final_err_offset, ...
            net_results(s).offset_bias, net_results(s).U_csac_err, rej_pct);
    end

    % Append this network's scenarios to long-form aggregate table.
    baseline = net_results(1).final_err_offset;
    pos_min = min(junction_positions);
    pos_range = max(junction_positions) - pos_min;
    if pos_range <= eps
        pos_range = 1;
    end

    for s = 1:num_scenarios
        r = net_results(s);
        sensor_idx = scenario_defs(s).sensor_indices;
        sensor_csac_ids = csac_ids(sensor_idx);
        sensor_pos_norm = (junction_positions(sensor_idx) - pos_min) / pos_range;

        one_sensor_zone = 'n/a';
        two_sensor_pattern = 'n/a';
        strategy_group = 'other';

        if isempty(sensor_idx)
            strategy_group = 'zero_sensors';
        elseif numel(sensor_idx) == 1
            one_sensor_zone = classify_single_sensor_zone(sensor_pos_norm(1));
            strategy_group = ['one_sensor_' one_sensor_zone];
        elseif numel(sensor_idx) == 2
            two_sensor_pattern = classify_two_sensor_pattern(sort(sensor_pos_norm(:)'));
            strategy_group = two_sensor_pattern;
        else
            strategy_group = sprintf('multi_sensor_%d', numel(sensor_idx));
        end

        improvement = 100 * (baseline - r.final_err_offset) / max(abs(baseline), eps);

        row = table();
        row.network_id = network_id;
        row.num_csacs = num_csacs;
        row.scenario = {r.name};
        row.num_sensors = numel(sensor_idx);
        row.sensor_csac_ids = {strjoin(arrayfun(@(x) sprintf('%d', x), sensor_csac_ids, 'UniformOutput', false), ',')};
        row.strategy_group = {strategy_group};
        row.one_sensor_zone = {one_sensor_zone};
        row.two_sensor_pattern = {two_sensor_pattern};
        row.tw_mae_offset = r.tw_mae_offset;
        row.final_err_offset = r.final_err_offset;
        row.offset_bias = r.offset_bias;
        row.offset_std = r.offset_std;
        row.tw_mae_U = r.tw_mae_U;
        row.final_err_U = r.final_err_U;
        row.U_csac_err = r.U_csac_err;
        row.rej_pct = r.rej_pct;
        row.mean_per_csac_final_err = mean(r.per_csac_final_err, 'omitnan');
        row.mean_per_csac_bias = mean(r.per_csac_bias, 'omitnan');
        row.improvement_vs_no_sensor_pct = improvement;

        if isempty(all_scenario_rows)
            all_scenario_rows = row;
        else
            all_scenario_rows = [all_scenario_rows; row]; %#ok<AGROW>
        end
    end

    %% ============================================================
    %% SAVE PER-NETWORK RESULTS
    %% ============================================================
    % Build summary table for this network
    summary = table();
    for s = 1:num_scenarios
        r = net_results(s);
        row = table();
        row.scenario = {r.name};
        row.num_sensors = r.num_sensors;
        row.tw_mae_offset = r.tw_mae_offset;
        row.final_err_offset = r.final_err_offset;
        row.offset_bias = r.offset_bias;
        row.offset_std = r.offset_std;
        row.tw_mae_U = r.tw_mae_U;
        row.final_err_U = r.final_err_U;
        row.U_csac_err = r.U_csac_err;
        row.rej_pct = r.rej_pct;

        % Add per-CSAC columns dynamically
        for c = 1:num_csacs
            row.(sprintf('csac_%d_final_err', csac_ids(c))) = r.per_csac_final_err(c);
            row.(sprintf('csac_%d_bias', csac_ids(c))) = r.per_csac_bias(c);
        end

        summary = [summary; row]; %#ok<AGROW>
    end

    writetable(summary, fullfile(output_folder, sprintf('network_%d_results.csv', network_id)));

    % Also save terminal-style text output
    fid = fopen(fullfile(output_folder, sprintf('network_%d_terminal.txt', network_id)), 'w');

    fprintf(fid, 'NETWORK %d RESULTS\n', network_id);
    fprintf(fid, 'Number of CSACs: %d\n', num_csacs);
    fprintf(fid, 'CSAC IDs: [%s]\n', strjoin(arrayfun(@(x) sprintf('%d', x), csac_ids, 'UniformOutput', false), ', '));
    fprintf(fid, 'Number of houses: %d\n', sum(cellfun(@(cs) cs.num_houses, all_cs)));
    fprintf(fid, 'U_csac true: %.4f W/m/K\n', U_csac_true);
    fprintf(fid, 'U_main true: %.4f W/m/K\n', U_main_true);
    fprintf(fid, '\n');

    % Main comparison table
    fprintf(fid, '%-20s | %8s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
        'Scenario', 'Off MAE', 'Off Fin', 'Off Bias', 'Off Std', 'U MAE', 'Uc err', 'Rej %%');
    fprintf(fid, '%s\n', repmat('-', 1, 100));
    for s = 1:num_scenarios
        r = net_results(s);
        fprintf(fid, '%-20s | %8.3f | %8.3f | %+8.3f | %8.3f | %8.4f | %+8.4f | %7.1f\n', ...
            r.name, r.tw_mae_offset, r.final_err_offset, r.offset_bias, ...
            r.offset_std, r.tw_mae_U, r.U_csac_err, r.rej_pct);
    end
    fprintf(fid, '\n');

    % Per-CSAC final offset error
    fprintf(fid, 'PER-CSAC FINAL OFFSET ERROR\n');
    fprintf(fid, '%-20s', 'Scenario');
    for c = 1:num_csacs
        fprintf(fid, ' | CSAC%d', csac_ids(c));
    end
    fprintf(fid, '\n%s\n', repmat('-', 1, 20 + num_csacs * 10));
    for s = 1:num_scenarios
        r = net_results(s);
        fprintf(fid, '%-20s', r.name);
        for c = 1:num_csacs
            fprintf(fid, ' | %7.3f', r.per_csac_final_err(c));
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');

    % Per-CSAC bias
    fprintf(fid, 'PER-CSAC OFFSET BIAS (signed mean)\n');
    fprintf(fid, '%-20s', 'Scenario');
    for c = 1:num_csacs
        fprintf(fid, ' | CSAC%d', csac_ids(c));
    end
    fprintf(fid, '\n%s\n', repmat('-', 1, 20 + num_csacs * 10));
    for s = 1:num_scenarios
        r = net_results(s);
        fprintf(fid, '%-20s', r.name);
        for c = 1:num_csacs
            fprintf(fid, ' | %+7.3f', r.per_csac_bias(c));
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');

    % Improvement table
    fprintf(fid, 'IMPROVEMENT vs NO SENSOR (Final offset error)\n');
    fprintf(fid, '%-20s | %10s | %10s\n', 'Scenario', 'Final Err', 'Improvement');
    fprintf(fid, '%s\n', repmat('-', 1, 45));
    baseline = net_results(1).final_err_offset;
    for s = 1:num_scenarios
        fin = net_results(s).final_err_offset;
        improvement = 100 * (baseline - fin) / baseline;
        fprintf(fid, '%-20s | %10.3f | %+9.1f%%\n', net_results(s).name, fin, improvement);
    end
    fprintf(fid, '\n');

    fclose(fid);

    fprintf('\nNetwork %d results saved to:\n', network_id);
    fprintf('  %s\n', fullfile(output_folder, sprintf('network_%d_results.csv', network_id)));
    fprintf('  %s\n', fullfile(output_folder, sprintf('network_%d_terminal.txt', network_id)));

    % Print to terminal too
    fprintf('\n');
    fprintf('NETWORK %d SUMMARY\n', network_id);
    fprintf('%-20s | %8s | %8s | %10s\n', 'Scenario', 'Off Fin', 'Off MAE', 'Improvement');
    fprintf('%s\n', repmat('-', 1, 55));
    for s = 1:num_scenarios
        r = net_results(s);
        improvement = 100 * (baseline - r.final_err_offset) / baseline;
        fprintf('%-20s | %8.3f | %8.3f | %+9.1f%%\n', ...
            r.name, r.final_err_offset, r.tw_mae_offset, improvement);
    end
    fprintf('\n');
end

%% ============================================================
%% SAVE OVERALL (CROSS-NETWORK) CONCLUSIONS
%% ============================================================
if isempty(all_scenario_rows)
    fprintf('\nNo scenario rows were collected; skipping overall summaries.\n');
else
    writetable(all_scenario_rows, fullfile(output_folder, 'overall_scenario_runs.csv'));

    [g_count, sensor_count_values] = findgroups(all_scenario_rows.num_sensors);
    by_sensor_count = table();
    by_sensor_count.num_sensors = sensor_count_values;
    by_sensor_count.n_runs = splitapply(@numel, all_scenario_rows.final_err_offset, g_count);
    by_sensor_count.mean_final_err_offset = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.final_err_offset, g_count);
    by_sensor_count.std_final_err_offset = splitapply(@(x) std(x, 'omitnan'), all_scenario_rows.final_err_offset, g_count);
    by_sensor_count.mean_improvement_vs_no_sensor_pct = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_count);
    by_sensor_count.pct_runs_better_than_no_sensor = splitapply(@(x) 100 * mean(x > 0, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_count);

    [g_one, one_zone_values] = findgroups(all_scenario_rows.one_sensor_zone);
    by_one_sensor_zone = table();
    by_one_sensor_zone.one_sensor_zone = one_zone_values;
    by_one_sensor_zone.n_runs = splitapply(@numel, all_scenario_rows.final_err_offset, g_one);
    by_one_sensor_zone.mean_final_err_offset = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.final_err_offset, g_one);
    by_one_sensor_zone.mean_improvement_vs_no_sensor_pct = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_one);
    by_one_sensor_zone.pct_runs_better_than_no_sensor = splitapply(@(x) 100 * mean(x > 0, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_one);
    by_one_sensor_zone = by_one_sensor_zone(~strcmp(by_one_sensor_zone.one_sensor_zone, 'n/a'), :);

    [g_two, two_pattern_values] = findgroups(all_scenario_rows.two_sensor_pattern);
    by_two_sensor_pattern = table();
    by_two_sensor_pattern.two_sensor_pattern = two_pattern_values;
    by_two_sensor_pattern.n_runs = splitapply(@numel, all_scenario_rows.final_err_offset, g_two);
    by_two_sensor_pattern.mean_final_err_offset = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.final_err_offset, g_two);
    by_two_sensor_pattern.mean_improvement_vs_no_sensor_pct = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_two);
    by_two_sensor_pattern.pct_runs_better_than_no_sensor = splitapply(@(x) 100 * mean(x > 0, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_two);
    by_two_sensor_pattern = by_two_sensor_pattern(~strcmp(by_two_sensor_pattern.two_sensor_pattern, 'n/a'), :);

    [g_strategy, strategy_values] = findgroups(all_scenario_rows.strategy_group);
    by_strategy = table();
    by_strategy.strategy_group = strategy_values;
    by_strategy.n_runs = splitapply(@numel, all_scenario_rows.final_err_offset, g_strategy);
    by_strategy.mean_final_err_offset = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.final_err_offset, g_strategy);
    by_strategy.std_final_err_offset = splitapply(@(x) std(x, 'omitnan'), all_scenario_rows.final_err_offset, g_strategy);
    by_strategy.mean_improvement_vs_no_sensor_pct = splitapply(@(x) mean(x, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_strategy);
    by_strategy.pct_runs_better_than_no_sensor = splitapply(@(x) 100 * mean(x > 0, 'omitnan'), all_scenario_rows.improvement_vs_no_sensor_pct, g_strategy);
    by_strategy = sortrows(by_strategy, 'mean_final_err_offset', 'ascend');

    writetable(by_sensor_count, fullfile(output_folder, 'overall_by_sensor_count.csv'));
    writetable(by_one_sensor_zone, fullfile(output_folder, 'overall_one_sensor_placement.csv'));
    writetable(by_two_sensor_pattern, fullfile(output_folder, 'overall_two_sensor_placement.csv'));
    writetable(by_strategy, fullfile(output_folder, 'overall_results.csv'));

    fid = fopen(fullfile(output_folder, 'overall_terminal.txt'), 'w');
    fprintf(fid, 'OVERALL STRATEGY RESULTS ACROSS %d NETWORKS\n\n', num_networks);

    fprintf(fid, 'A) SENSOR COUNT (0 vs 1 vs 2)\n');
    fprintf(fid, '%-10s | %6s | %10s | %10s | %12s\n', ...
        'Sensors', 'Runs', 'OffFin(mu)', 'OffFin(sd)', 'BetterThan0');
    fprintf(fid, '%s\n', repmat('-', 1, 62));
    for i = 1:height(by_sensor_count)
        fprintf(fid, '%-10d | %6d | %10.3f | %10.3f | %10.1f%%\n', ...
            by_sensor_count.num_sensors(i), by_sensor_count.n_runs(i), ...
            by_sensor_count.mean_final_err_offset(i), by_sensor_count.std_final_err_offset(i), ...
            by_sensor_count.pct_runs_better_than_no_sensor(i));
    end
    fprintf(fid, '\n');

    fprintf(fid, 'B) ONE SENSOR PLACEMENT (beginning / middle / end)\n');
    fprintf(fid, '%-18s | %6s | %10s | %12s\n', ...
        'Zone', 'Runs', 'OffFin(mu)', 'BetterThan0');
    fprintf(fid, '%s\n', repmat('-', 1, 58));
    for i = 1:height(by_one_sensor_zone)
        fprintf(fid, '%-18s | %6d | %10.3f | %10.1f%%\n', ...
            by_one_sensor_zone.one_sensor_zone{i}, by_one_sensor_zone.n_runs(i), ...
            by_one_sensor_zone.mean_final_err_offset(i), ...
            by_one_sensor_zone.pct_runs_better_than_no_sensor(i));
    end
    fprintf(fid, '\n');

    fprintf(fid, 'C) TWO SENSOR PLACEMENT (far / close-begin / close-end)\n');
    fprintf(fid, '%-22s | %6s | %10s | %12s\n', ...
        'Pattern', 'Runs', 'OffFin(mu)', 'BetterThan0');
    fprintf(fid, '%s\n', repmat('-', 1, 64));
    for i = 1:height(by_two_sensor_pattern)
        fprintf(fid, '%-22s | %6d | %10.3f | %10.1f%%\n', ...
            by_two_sensor_pattern.two_sensor_pattern{i}, by_two_sensor_pattern.n_runs(i), ...
            by_two_sensor_pattern.mean_final_err_offset(i), ...
            by_two_sensor_pattern.pct_runs_better_than_no_sensor(i));
    end
    fprintf(fid, '\n');

    fprintf(fid, 'D) STRATEGY RANKING\n');
    fprintf(fid, '%-24s | %6s | %10s | %10s | %12s\n', ...
        'Strategy', 'Runs', 'OffFin(mu)', 'OffFin(sd)', 'BetterThan0');
    fprintf(fid, '%s\n', repmat('-', 1, 74));
    for i = 1:height(by_strategy)
        fprintf(fid, '%-24s | %6d | %10.3f | %10.3f | %10.1f%%\n', ...
            by_strategy.strategy_group{i}, by_strategy.n_runs(i), ...
            by_strategy.mean_final_err_offset(i), by_strategy.std_final_err_offset(i), ...
            by_strategy.pct_runs_better_than_no_sensor(i));
    end
    fclose(fid);

    fprintf('Overall results saved to:\n');
    fprintf('  %s\n', fullfile(output_folder, 'overall_scenario_runs.csv'));
    fprintf('  %s\n', fullfile(output_folder, 'overall_by_sensor_count.csv'));
    fprintf('  %s\n', fullfile(output_folder, 'overall_one_sensor_placement.csv'));
    fprintf('  %s\n', fullfile(output_folder, 'overall_two_sensor_placement.csv'));
    fprintf('  %s\n', fullfile(output_folder, 'overall_results.csv'));
    fprintf('  %s\n', fullfile(output_folder, 'overall_terminal.txt'));
end

fprintf('\nAll networks processed. Results saved to results/sensor_test/\n');

function zone = classify_single_sensor_zone(pos_norm)
if pos_norm <= 1/3
    zone = 'beginning';
elseif pos_norm >= 2/3
    zone = 'end';
else
    zone = 'middle';
end
end

function pattern = classify_two_sensor_pattern(pos_norm_sorted)
spacing = abs(pos_norm_sorted(2) - pos_norm_sorted(1));
pair_center = mean(pos_norm_sorted);

if spacing > 1/3
    pattern = 'two_far_apart';
else
    if pair_center <= 1/3
        pattern = 'two_close_beginning';
    elseif pair_center >= 2/3
        pattern = 'two_close_end';
    else
        pattern = 'two_close_middle';
    end
end
end