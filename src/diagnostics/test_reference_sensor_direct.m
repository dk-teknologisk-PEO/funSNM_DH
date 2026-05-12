% test_reference_sensor_direct.m
% Diagnostic variant of test_reference_sensor.m.
%
% Goal: isolate whether the main-pipe temperature-profile propagation in
% the original test is the source of the performance regression on
% network_2 - network_4 when reference sensors are installed.
%
% Two scenarios only:
%   1) no_sensor         - baseline (each CSAC on its own)
%   2) reference_direct  - reference sensor at EVERY main-pipe/CSAC junction.
%                          Reference temperature is taken NOISE-FREE from
%                          the data file (J_Main_<id>_s) and assigned
%                          directly to cs.T_inlet_from_main. T_inlet_fitted
%                          is completely bypassed; no propagation along the
%                          main pipe, no inverse-variance fusion.
%                          get_supply_temp.m still handles the
%                          main-junction -> service-pipe step inside
%                          process_csac_timestep.

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

networks_to_test = config.project.datasets.datasets(:)';
num_networks = numel(networks_to_test);

%% Scenario definitions
scenario_defs = struct();
scenario_defs(1).name = 'no_sensor';
scenario_defs(1).use_reference = false;
scenario_defs(2).name = 'reference_direct';
scenario_defs(2).use_reference = true;
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

    % Load NOISE-FREE reference temperature column for every CSAC junction.
    % These are the "ideal sensor" readings used in scenario 2.
    ref_data = cell(num_csacs, 1);
    for c = 1:num_csacs
        col_name = sprintf('J_Main_%d_s', csac_ids(c));
        if ismember(col_name, network_data.Properties.VariableNames)
            ref_data{c} = struct('timestamps', network_data.timestamp, ...
                'temperatures', network_data.(col_name), 'valid', true);
        else
            ref_data{c} = struct('timestamps', [], 'temperatures', [], 'valid', false);
            warning('Network %d: missing column %s', network_id, col_name);
        end
    end

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

        U_csac_cfg = config.project.csac_U_estimation;
        U_csac_update_counter = 0;

        use_reference = scenario_defs(s).use_reference;

        % Diagnostic log: T_inlet_fitted (meter-based) vs reference T.
        % Only meaningful when use_reference == true (we need both side-by-side).
        diag_log = struct('time', {}, 'csac_id', {}, 'T_fitted', {}, ...
            'T_sigma', {}, 'T_ref', {});

        % Timestep loop
        for t = 1:length(timestamps)
            time = timestamps(t);
            current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
            if isempty(current_T_soil_C)
                continue
            end

            % Step 1: Process each CSAC independently (always run -- this
            % gives the season_state and the baseline T_inlet_fitted).
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

            % Step 2: If reference sensors are used, OVERRIDE T_inlet_from_main
            % with the noise-free reference temperature at this CSAC's
            % junction (bypassing T_inlet_fitted entirely) and re-process.
            % BEFORE overriding, log the meter-based T_inlet_fitted vs the
            % reference temperature for diagnostic comparison.
            if use_reference
                for c = 1:num_csacs
                    if ~ref_data{c}.valid
                        continue
                    end
                    ref_idx = find(ref_data{c}.timestamps == time, 1);
                    if isempty(ref_idx)
                        continue
                    end
                    ref_T = ref_data{c}.temperatures(ref_idx);
                    if ~isfinite(ref_T)
                        continue
                    end

                    % --- Diagnostic: capture meter-based estimate from step 1 ---
                    cs_c = all_cs{c};
                    if isfield(cs_c, 'T_inlet_fitted') && isfinite(cs_c.T_inlet_fitted)
                        entry = struct( ...
                            'time',    time, ...
                            'csac_id', csac_ids(c), ...
                            'T_fitted', cs_c.T_inlet_fitted, ...
                            'T_sigma',  cs_c.T_inlet_sigma, ...
                            'T_ref',    ref_T);
                        diag_log(end+1) = entry; %#ok<AGROW>
                    end

                    all_cs{c}.T_inlet_from_main = ref_T;
                end

                % Re-process CSACs with the forced T_inlet_from_main
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
                fullfile('results', 'sensor_test_direct', sprintf('net%d', network_id), ...
                scenario_defs(s).name), kpi_config, all_true_traj{c}, network_id);
            all_kpi = [all_kpi; kpi]; %#ok<AGROW>
        end

        % Additional metrics
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

        max_csacs = 4;
        padded_err = nan(1, max_csacs);
        padded_bias = nan(1, max_csacs);
        padded_err(1:num_csacs) = per_csac_final_err;
        padded_bias(1:num_csacs) = per_csac_bias;
        agg_results(s).all_per_csac_final_err = [agg_results(s).all_per_csac_final_err; padded_err];
        agg_results(s).all_per_csac_bias = [agg_results(s).all_per_csac_bias; padded_bias];

        fprintf('  Off MAE=%.3f, Final=%.3f, Bias=%+.3f, Uc_err=%+.4f, Rej=%.1f%%\n', ...
            mean(all_kpi.tw_mae_offset, 'omitnan'), ...
            mean(abs(all_kpi.final_err_offset), 'omitnan'), ...
            mean(offset_errors), ...
            shared_U_csac - U_csac_true, rej_pct);

        %% Diagnostic: save & summarize T_inlet_fitted vs reference T
        if use_reference && ~isempty(diag_log)
            diag_tbl = struct2table(diag_log);
            diag_tbl.diff = diag_tbl.T_fitted - diag_tbl.T_ref;

            diag_folder = fullfile('results', 'sensor_test_direct', ...
                sprintf('net%d', network_id), scenario_defs(s).name);
            if ~exist(diag_folder, 'dir'), mkdir(diag_folder); end
            writetable(diag_tbl, fullfile(diag_folder, 'T_inlet_fitted_vs_ref.csv'));

            % Summary: overall and tail (last 30%% of samples, after warm-up)
            fprintf('\n  T_inlet_fitted vs T_ref (per CSAC):\n');
            fprintf('  %-7s | %6s | %8s | %8s | %8s | %8s\n', ...
                'CSAC', 'N', 'meanDiff', 'stdDiff', 'tail_md', 'tail_sd');
            for c = 1:num_csacs
                mask_all = diag_tbl.csac_id == csac_ids(c);
                d_all = diag_tbl.diff(mask_all);
                if isempty(d_all), continue; end
                n = numel(d_all);
                tail_start = max(1, round(0.7 * n));
                d_tail = d_all(tail_start:end);
                fprintf('  %-7d | %6d | %+8.3f | %8.3f | %+8.3f | %8.3f\n', ...
                    csac_ids(c), n, mean(d_all), std(d_all), ...
                    mean(d_tail), std(d_tail));
            end

            % --- Plot: T_inlet_fitted - T_ref vs time, one line per CSAC ---
            fig = figure('Visible', 'off', 'Position', [100 100 1000 500]);
            hold on; grid on;
            legend_entries = {};
            for c = 1:num_csacs
                mask_c = diag_tbl.csac_id == csac_ids(c);
                if ~any(mask_c), continue; end
                plot(diag_tbl.time(mask_c), diag_tbl.diff(mask_c), '-', 'LineWidth', 1);
                legend_entries{end+1} = sprintf('CSAC %d', csac_ids(c)); %#ok<AGROW>
            end
            yline(0, 'k--', 'HandleVisibility', 'off');
            xlabel('Time');
            ylabel('T_{fitted} - T_{ref} [\circC]');
            title(sprintf('Network %d: meter-based inlet temp vs reference', network_id));
            legend(legend_entries, 'Location', 'eastoutside');
            hold off;

            fig_path = fullfile(diag_folder, ...
                sprintf('T_inlet_fitted_vs_ref_error_net%d', network_id));
            save_figure(fig, fig_path);
            close(fig);
        end
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

%% Per-CSAC breakdown
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

%% Improvement summary vs baseline
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
output_folder = fullfile('results', 'sensor_test_direct');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

summary = table();
for s = 1:num_scenarios
    a = agg_results(s);
    row = table();
    row.scenario = {a.name};
    row.use_reference = scenario_defs(s).use_reference;
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
writetable(summary, fullfile(output_folder, 'sensor_comparison_direct.csv'));

fprintf('\nResults saved to %s\n', fullfile(output_folder, 'sensor_comparison_direct.csv'));
fprintf('Test complete.\n');
