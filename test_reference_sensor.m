% test_reference_sensor.m
% Tests the impact of reference sensor placement on estimation accuracy.

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

%% Load data
networks = config.project.datasets.datasets(:)';
network_id = networks(1);
[meter_data, network_data, topology] = importData(config, network_id);
timestamps = unique(meter_data.timestamp);
csac_ids = [topology.cul_de_sacs.id];
num_csacs = numel(csac_ids);

%% True U values for reference
U_csac_true = topology.pipe_parameters.csac_pipe.insulation_W_m_K;
U_main_true = topology.pipe_parameters.main_pipe.insulation_W_m_K;

%% Define sensor scenarios
sensor_noise = 0.1; % °C

scenarios = struct();

scenarios(1).name = 'no_sensor';
scenarios(1).sensor = struct('type', 'none', 'csac_id', 0, 'noise_std', 0, 'seed', 42);
scenarios(1).coupling_enabled = false;

scenarios(2).name = 'plant_inlet';
scenarios(2).sensor = struct('type', 'plant', 'csac_id', 0, 'noise_std', sensor_noise, 'seed', 42);
scenarios(2).coupling_enabled = true;

for c = 1:num_csacs
    idx = 2 + c;
    scenarios(idx).name = sprintf('csac_%d_junction', csac_ids(c));
    scenarios(idx).sensor = struct('type', 'csac_junction', 'csac_id', csac_ids(c), ...
        'noise_std', sensor_noise, 'seed', 42);
    scenarios(idx).coupling_enabled = true;
end

%% Run each scenario
results = struct();

for s = 1:numel(scenarios)
    fprintf('\n========================================\n');
    fprintf('SCENARIO %d/%d: %s\n', s, numel(scenarios), scenarios(s).name);
    fprintf('========================================\n');

    % Initialize all CSACs
    [all_cs, all_true_traj] = initialize_all_csacs(topology, meter_data, timestamps, ...
        Q_base, R_base, P_base, innovation_gate_initial, config, network_id, drift_config);

    % Initialize shared U_csac (same random init for all scenarios)
    shared_U_csac = all_cs{1}.U_csac;
    for c = 1:num_csacs
        all_cs{c}.U_csac = shared_U_csac;
    end

    % U_main: use true value with small uncertainty for propagation
    % NOT randomly initialized — we want to isolate the sensor effect
    shared_U_main = U_main_true;
    U_main_uncertainty = 0.01; % std for propagation uncertainty

    U_csac_cfg = config.project.csac_U_estimation;
    U_main_cfg = config.project.main_pipe_U_estimation;
    U_csac_update_counter = 0;
    U_csac_history = shared_U_csac;

    % Load reference sensor data
    ref_data = load_reference_sensor_data(network_data, scenarios(s).sensor, network_id);

    % Timestep loop
    main_coupling_counter = 0;

    for t = 1:length(timestamps)
        time = timestamps(t);
        current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
        if isempty(current_T_soil_C)
            continue
        end

        % Get reference temperature for this timestep
        ref_T = NaN;
        if ref_data.valid
            ref_idx = find(ref_data.timestamps == time, 1);
            if ~isempty(ref_idx)
                ref_T = ref_data.temperatures(ref_idx);
            end
        end

        % Process each CSAC
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

        % Main pipe coupling using sensor data with uncertainty propagation
        main_coupling_counter = main_coupling_counter + 1;
        if scenarios(s).coupling_enabled && isfinite(ref_T) && ...
           main_coupling_counter >= config.project.main_pipe_coupling.warmup_timesteps

            % Get junction positions and flows
            junction_positions = zeros(num_csacs, 1);
            junction_flows = zeros(num_csacs, 1);
            for c = 1:num_csacs
                topo_idx = find([topology.cul_de_sacs.id] == csac_ids(c));
                junction_positions(c) = topology.cul_de_sacs(topo_idx).dist_on_main_m;
                junction_flows(c) = all_cs{c}.current_total_flow;
            end

            % Determine sensor position
            if strcmp(scenarios(s).sensor.type, 'plant')
                sensor_pos = 0;
            else
                topo_idx = find([topology.cul_de_sacs.id] == scenarios(s).sensor.csac_id);
                sensor_pos = topology.cul_de_sacs(topo_idx).dist_on_main_m;
            end

            % Propagate sensor temperature with uncertainty
            [T_prop, T_unc] = propagate_sensor_temperature(...
                ref_T, sensor_pos, scenarios(s).sensor.noise_std, ...
                junction_positions, junction_flows, ...
                shared_U_main, U_main_uncertainty, current_T_soil_C, topology);

            % Blend propagated temperature with CSAC's own estimate
            for c = 1:num_csacs
                if isfinite(T_prop(c)) && isfinite(T_unc(c)) && T_unc(c) > 0
                    T_sensor = T_prop(c);
                    sigma_sensor = T_unc(c);

                    if isfield(all_cs{c}, 'T_inlet_fitted') && isfinite(all_cs{c}.T_inlet_fitted)
                        T_own = all_cs{c}.T_inlet_fitted;
                        sigma_own = 0.5;

                        w_sensor = 1 / sigma_sensor^2;
                        w_own = 1 / sigma_own^2;
                        T_blended = (w_sensor * T_sensor + w_own * T_own) / (w_sensor + w_own);

                        all_cs{c}.T_inlet_from_main = T_blended;
                    end
                end
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
                U_csac_history(end+1) = shared_U_csac;
            end
        end
    end

    % Compute KPIs
    output_folder = fullfile('results', 'sensor_test', scenarios(s).name);
    all_kpi = table();
    for c = 1:num_csacs
        kpi = compute_and_save_network_kpis(all_cs{c}, csac_ids(c), output_folder, ...
            kpi_config, all_true_traj{c}, network_id);
        all_kpi = [all_kpi; kpi]; %#ok<AGROW>
    end

    % Compute additional metrics
    % Mean offset bias (absolute level shift across all houses)
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
    mean_offset_bias = mean(offset_errors);
    std_offset_error = std(offset_errors);

    results(s).name = scenarios(s).name;
    results(s).kpis = all_kpi;
    results(s).U_csac_final = shared_U_csac;
    results(s).U_csac_history = U_csac_history;
    results(s).mean_offset_bias = mean_offset_bias;
    results(s).std_offset_error = std_offset_error;
    results(s).all_cs = all_cs;

    % Find U_csac convergence step (first time within 0.005 of true)
    U_csac_converged_step = NaN;
    for k = 1:numel(U_csac_history)
        if abs(U_csac_history(k) - U_csac_true) < 0.005
            U_csac_converged_step = k;
            break;
        end
    end
    results(s).U_csac_converged_step = U_csac_converged_step;

    fprintf('  Final: U_csac=%.4f (true=%.4f, err=%.4f)\n', ...
        shared_U_csac, U_csac_true, shared_U_csac - U_csac_true);
    fprintf('  Offset bias: mean=%.4f°C, std=%.4f°C\n', mean_offset_bias, std_offset_error);
    fprintf('  U_csac converged at step: %s\n', mat2str(U_csac_converged_step));
end

%% Print comparison table
fprintf('\n========================================\n');
fprintf('REFERENCE SENSOR COMPARISON\n');
fprintf('========================================\n');
fprintf('%-20s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'Scenario', 'Off MAE', 'Off Fin', 'Off Bias', 'Off Std', 'U MAE', 'Uc err', 'Uc conv', 'Rej %%');
fprintf('%s\n', repmat('-', 1, 105));

for s = 1:numel(results)
    k = results(s).kpis;
    total_accept = sum(cellfun(@(cs) cs.gate_accept_count, results(s).all_cs));
    total_reject = sum(cellfun(@(cs) cs.gate_reject_count, results(s).all_cs));
    rej_pct = 100 * total_reject / max(1, total_accept + total_reject);

    conv_str = sprintf('%d', results(s).U_csac_converged_step);
    if isnan(results(s).U_csac_converged_step)
        conv_str = 'N/A';
    end

    fprintf('%-20s | %8.3f | %8.3f | %+8.3f | %8.3f | %8.4f | %+8.4f | %8s | %7.1f\n', ...
        results(s).name, ...
        mean(k.tw_mae_offset, 'omitnan'), ...
        mean(abs(k.final_err_offset), 'omitnan'), ...
        results(s).mean_offset_bias, ...
        results(s).std_offset_error, ...
        mean(k.tw_mae_U, 'omitnan'), ...
        results(s).U_csac_final - U_csac_true, ...
        conv_str, ...
        rej_pct);
end
fprintf('========================================\n');

%% Per-CSAC breakdown
fprintf('\nPER-CSAC FINAL OFFSET ERROR (absolute)\n');
fprintf('%-20s', 'Scenario');
for c = 1:num_csacs
    fprintf(' | CSAC%d', csac_ids(c));
end
fprintf('\n%s\n', repmat('-', 1, 20 + num_csacs * 10));

for s = 1:numel(results)
    fprintf('%-20s', results(s).name);
    k = results(s).kpis;
    for c = 1:num_csacs
        csac_mask = k.csac_id == csac_ids(c);
        fprintf(' | %7.3f', mean(abs(k.final_err_offset(csac_mask)), 'omitnan'));
    end
    fprintf('\n');
end
fprintf('\n');

%% Per-CSAC offset bias (signed mean error)
fprintf('PER-CSAC OFFSET BIAS (signed mean)\n');
fprintf('%-20s', 'Scenario');
for c = 1:num_csacs
    fprintf(' | CSAC%d', csac_ids(c));
end
fprintf('\n%s\n', repmat('-', 1, 20 + num_csacs * 10));

for s = 1:numel(results)
    fprintf('%-20s', results(s).name);
    for c = 1:num_csacs
        cs = results(s).all_cs{c};
        offsets = zeros(cs.num_houses, 1);
        true_off = cs.ground_truth.true_offset;
        for i = 1:cs.num_houses
            offsets(i) = cs.ukf_states{i}.x(1);
        end
        bias = mean(offsets - true_off);
        fprintf(' | %+7.3f', bias);
    end
    fprintf('\n');
end
fprintf('========================================\n');

%% Convergence speed comparison
fprintf('\nCONVERGENCE SPEED\n');
fprintf('%-20s | %12s | %12s\n', 'Scenario', 'Uc conv step', 'Off conv days');
fprintf('%s\n', repmat('-', 1, 50));
for s = 1:numel(results)
    k = results(s).kpis;
    mean_conv = mean(k.conv_days_offset, 'omitnan');
    conv_str = sprintf('%d', results(s).U_csac_converged_step);
    if isnan(results(s).U_csac_converged_step)
        conv_str = 'N/A';
    end
    fprintf('%-20s | %12s | %12.0f\n', results(s).name, conv_str, mean_conv);
end
fprintf('========================================\n');

%% Save results
output_folder = fullfile('results', 'sensor_test');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

comparison = table();
for s = 1:numel(results)
    k = results(s).kpis;
    row = table();
    row.scenario = {results(s).name};
    row.mean_tw_mae_offset = mean(k.tw_mae_offset, 'omitnan');
    row.mean_final_err_offset = mean(abs(k.final_err_offset), 'omitnan');
    row.offset_bias = results(s).mean_offset_bias;
    row.offset_std = results(s).std_offset_error;
    row.mean_tw_mae_U = mean(k.tw_mae_U, 'omitnan');
    row.U_csac_error = results(s).U_csac_final - U_csac_true;
    row.U_csac_conv_step = results(s).U_csac_converged_step;
    comparison = [comparison; row]; %#ok<AGROW>
end
writetable(comparison, fullfile(output_folder, 'sensor_comparison.csv'));

fprintf('\nSensor test complete. Results saved to results/sensor_test/\n');