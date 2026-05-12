% test_no_sensor_inlet_error.m
%
% Quantify the estimation error of the CSAC-pipe inlet temperature in the
% no_sensor scenario (each CSAC runs on its own, no reference sensor on the
% main pipe). For every CSAC at every timestep we log:
%
%   T_fitted = cs.T_inlet_fitted   (meter-based estimate of the CSAC-pipe
%                                   inlet temperature)
%   T_ref    = J_Main_<id>_s       (noise-free reference temperature from
%                                   the simulation data file)
%   err      = T_fitted - T_ref
%
% All samples are stored in a single timetable. After the main loop the
% data are split into heating seasons (1 July -> 1 July) and a histogram +
% mean/std is produced for each season.
%
% Derived from test_reference_sensor_direct.m -- the two reference-sensor
% scenarios have been stripped out; only the no_sensor baseline remains.

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

% Pre-scan topologies for max CSAC count (used when padding result rows).
max_csacs = 0;
topology_path = config.project.paths.topology;
for pre_idx = 1:num_networks
    network_ID_str = sprintf('%02d', networks_to_test(pre_idx));
    topo_file = fullfile(topology_path, strcat("network_", network_ID_str, ".json"));
    pre_topology = jsondecode(fileread(topo_file));
    max_csacs = max(max_csacs, numel(pre_topology.cul_de_sacs));
end
fprintf('Maximum CSACs across networks: %d\n', max_csacs);

fprintf('Running no_sensor scenario across %d networks\n', num_networks);

% Cross-network aggregator for the inlet-temperature error.
% Each row: timestamp, network_id, csac_id, T_fitted, T_ref, err
no_sensor_inlet_log = table('Size', [0 6], ...
    'VariableTypes', {'datetime','double','double','double','double','double'}, ...
    'VariableNames', {'time','network_id','csac_id','T_fitted','T_ref','err'});

% Per-network aggregate KPIs (single scenario -> flat structure).
agg = struct( ...
    'all_tw_mae_offset',     [], ...
    'all_final_err_offset',  [], ...
    'all_tw_mae_U',          [], ...
    'all_final_err_U',       [], ...
    'all_offset_bias',       [], ...
    'all_offset_std',        [], ...
    'all_U_csac_err',        [], ...
    'all_rej_pct',           [], ...
    'all_per_csac_final_err',[], ...
    'all_per_csac_bias',     []);

%% Main loop over networks
for net_idx = 1:num_networks
    network_id = networks_to_test(net_idx);
    fprintf('\n############################################################\n');
    fprintf('NETWORK %d (%d/%d) - no_sensor\n', network_id, net_idx, num_networks);
    fprintf('############################################################\n');

    [meter_data, network_data, topology] = importData(config, network_id);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];
    num_csacs = numel(csac_ids);

    U_csac_true = topology.pipe_parameters.csac_pipe.insulation_W_m_K;

    % Load NOISE-FREE reference temperature column for every CSAC junction.
    % Used only for error analysis -- never fed into the filter.
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

    % Initialize all CSACs
    [all_cs, all_true_traj] = initialize_all_csacs(topology, meter_data, timestamps, ...
        Q_base, R_base, P_base, innovation_gate_initial, config, network_id, drift_config);

    shared_U_csac = all_cs{1}.U_csac;
    for c = 1:num_csacs
        all_cs{c}.U_csac = shared_U_csac;
    end

    U_csac_cfg = config.project.csac_U_estimation;
    U_csac_update_counter = 0;

    % Per-timestep log of T_inlet_fitted vs noise-free reference T.
    diag_log = struct('time', {}, 'csac_id', {}, 'T_fitted', {}, ...
        'T_sigma', {}, 'T_ref', {});

    % Timestep loop
    for t = 1:length(timestamps)
        time = timestamps(t);
        current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
        if isempty(current_T_soil_C)
            continue
        end

        % Step 1: Process each CSAC independently.
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

        % Step 1b: log T_inlet_fitted vs reference T (one entry per CSAC
        % per timestep). Reference column is the noise-free J_Main_<id>_s
        % from the simulation data file.
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
            cs_c = all_cs{c};
            if isfield(cs_c, 'T_inlet_fitted') && isfinite(cs_c.T_inlet_fitted)
                entry = struct( ...
                    'time',     time, ...
                    'csac_id',  csac_ids(c), ...
                    'T_fitted', cs_c.T_inlet_fitted, ...
                    'T_sigma',  cs_c.T_inlet_sigma, ...
                    'T_ref',    ref_T);
                diag_log(end+1) = entry; %#ok<AGROW>
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
            fullfile('results', 'no_sensor_inlet_error', sprintf('net%d', network_id)), ...
            kpi_config, all_true_traj{c}, network_id);
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

    agg.all_tw_mae_offset(end+1)    = mean(all_kpi.tw_mae_offset, 'omitnan');
    agg.all_final_err_offset(end+1) = mean(abs(all_kpi.final_err_offset), 'omitnan');
    agg.all_tw_mae_U(end+1)         = mean(all_kpi.tw_mae_U, 'omitnan');
    agg.all_final_err_U(end+1)      = mean(abs(all_kpi.final_err_U), 'omitnan');
    agg.all_offset_bias(end+1)      = mean(offset_errors);
    agg.all_offset_std(end+1)       = std(offset_errors);
    agg.all_U_csac_err(end+1)       = shared_U_csac - U_csac_true;
    agg.all_rej_pct(end+1)          = rej_pct;

    padded_err = nan(1, max_csacs);
    padded_bias = nan(1, max_csacs);
    padded_err(1:num_csacs)  = per_csac_final_err;
    padded_bias(1:num_csacs) = per_csac_bias;
    agg.all_per_csac_final_err = [agg.all_per_csac_final_err; padded_err];
    agg.all_per_csac_bias      = [agg.all_per_csac_bias; padded_bias];

    fprintf('  Off MAE=%.3f, Final=%.3f, Bias=%+.3f, Uc_err=%+.4f, Rej=%.1f%%\n', ...
        mean(all_kpi.tw_mae_offset, 'omitnan'), ...
        mean(abs(all_kpi.final_err_offset), 'omitnan'), ...
        mean(offset_errors), ...
        shared_U_csac - U_csac_true, rej_pct);

    %% Save & summarize T_inlet_fitted vs reference T (this network)
    if ~isempty(diag_log)
        diag_tbl = struct2table(diag_log);
        diag_tbl.diff = diag_tbl.T_fitted - diag_tbl.T_ref;

        % Append to cross-network aggregator.
        add_tbl = table(diag_tbl.time, ...
            repmat(network_id, height(diag_tbl), 1), ...
            diag_tbl.csac_id, diag_tbl.T_fitted, diag_tbl.T_ref, ...
            diag_tbl.diff, ...
            'VariableNames', {'time','network_id','csac_id','T_fitted','T_ref','err'});
        no_sensor_inlet_log = [no_sensor_inlet_log; add_tbl]; %#ok<AGROW>

        diag_folder = fullfile('results', 'no_sensor_inlet_error', ...
            sprintf('net%d', network_id));
        if ~exist(diag_folder, 'dir'), mkdir(diag_folder); end
        writetable(diag_tbl, fullfile(diag_folder, 'T_inlet_fitted_vs_ref.csv'));

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

    %% Per-CSAC time-series plots: offset error and U_service error
    ts_folder = fullfile('results', 'no_sensor_inlet_error', sprintf('net%d', network_id));
    if ~exist(ts_folder, 'dir'), mkdir(ts_folder); end

    for c = 1:num_csacs
        cs = all_cs{c};
        log = cs.logger;
        t_axis = log.timestamps;
        valid_t = ~isnat(t_axis);
        if ~any(valid_t), continue; end

        true_off = cs.ground_truth.true_offset;
        true_U   = cs.ground_truth.true_U;

        fig = figure('Visible', 'off', 'Position', [100 100 1100 700]);

        subplot(2, 1, 1); hold on; grid on;
        leg = cell(1, cs.num_houses);
        for i = 1:cs.num_houses
            off_traj = squeeze(log.state_estimates(1, i, :));
            err = off_traj(:) - true_off(i);
            plot(t_axis(valid_t), err(valid_t), '-', 'LineWidth', 1);
            leg{i} = sprintf('H%d', cs.house_ids(i));
        end
        yline(0, 'k--', 'HandleVisibility', 'off');
        ylabel('Offset error [\circC]');
        title(sprintf('Net %d, CSAC %d, no\\_sensor — meter offset error', ...
            network_id, csac_ids(c)));
        legend(leg, 'Location', 'eastoutside');
        hold off;

        subplot(2, 1, 2); hold on; grid on;
        for i = 1:cs.num_houses
            U_traj = squeeze(log.state_estimates(2, i, :));
            err = U_traj(:) - true_U(i);
            plot(t_axis(valid_t), err(valid_t), '-', 'LineWidth', 1);
        end
        yline(0, 'k--', 'HandleVisibility', 'off');
        xlabel('Time');
        ylabel('U_{service} error [W/m/K]');
        title('Service-pipe U-value error');
        legend(leg, 'Location', 'eastoutside');
        hold off;

        fig_path = fullfile(ts_folder, ...
            sprintf('timeseries_net%d_csac%d', network_id, csac_ids(c)));
        save_figure(fig, fig_path);
        close(fig);
    end
end

%% Print aggregate KPI summary
fprintf('\n\n############################################################\n');
fprintf('AGGREGATE no_sensor RESULTS ACROSS %d NETWORKS\n', num_networks);
fprintf('############################################################\n\n');

fprintf('%-8s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'Off MAE', 'Off Fin', 'Off Bias', 'Off Std', 'U MAE', 'Uc err', 'Rej %%');
fprintf('%s\n', repmat('-', 1, 80));
fprintf('%8.3f | %8.3f | %+8.3f | %8.3f | %8.4f | %+8.4f | %7.1f\n', ...
    mean(agg.all_tw_mae_offset), ...
    mean(agg.all_final_err_offset), ...
    mean(agg.all_offset_bias), ...
    mean(agg.all_offset_std), ...
    mean(agg.all_tw_mae_U), ...
    mean(agg.all_U_csac_err), ...
    mean(agg.all_rej_pct));

fprintf('\nPER-CSAC FINAL OFFSET ERROR (averaged across %d networks)\n', num_networks);
avg_per_csac = mean(agg.all_per_csac_final_err, 1);
for k = 1:max_csacs
    fprintf('  CSAC %d: %8.3f\n', k-1, avg_per_csac(k));
end

%% Save aggregate KPI CSV
output_folder = fullfile('results', 'no_sensor_inlet_error');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

summary = table( ...
    {'no_sensor'}, ...
    mean(agg.all_tw_mae_offset), ...
    mean(agg.all_final_err_offset), ...
    mean(agg.all_offset_bias), ...
    mean(agg.all_offset_std), ...
    mean(agg.all_tw_mae_U), ...
    mean(agg.all_U_csac_err), ...
    mean(agg.all_rej_pct), ...
    num_networks, ...
    'VariableNames', {'scenario','mean_tw_mae_offset','mean_final_err_offset', ...
        'mean_offset_bias','mean_offset_std','mean_tw_mae_U', ...
        'mean_U_csac_err','mean_rej_pct','num_networks'});
writetable(summary, fullfile(output_folder, 'no_sensor_summary.csv'));
fprintf('\nKPI summary saved to %s\n', fullfile(output_folder, 'no_sensor_summary.csv'));

%% =====================================================================
%% Inlet-temperature error per heating season (1 July -> 1 July)
%% =====================================================================
if ~isempty(no_sensor_inlet_log)
    no_sensor_inlet_log = sortrows(no_sensor_inlet_log, 'time');
    inlet_err_tt = table2timetable(no_sensor_inlet_log, 'RowTimes', 'time');
    save(fullfile(output_folder, 'inlet_err_timetable.mat'), 'inlet_err_tt');
    writetimetable(inlet_err_tt, fullfile(output_folder, 'inlet_err_timetable.csv'));
    fprintf('\nInlet-temp error timetable: %d rows saved to %s\n', ...
        height(inlet_err_tt), fullfile(output_folder, 'inlet_err_timetable.csv'));

    t_all = inlet_err_tt.time;
    yr_min = year(min(t_all));
    yr_max = year(max(t_all));
    season_starts = datetime(yr_min-1:yr_max+1, 7, 1, 'TimeZone', t_all.TimeZone);

    fprintf('\nHeating-season inlet-temperature error (T_fitted - T_ref):\n');
    fprintf('%-25s | %8s | %+9s | %8s\n', 'Season (Jul-Jun)', 'N', 'mean [C]', 'std [C]');
    fprintf('%s\n', repmat('-', 1, 60));

    season_summary = table();
    for k = 1:numel(season_starts)-1
        s_start = season_starts(k);
        s_end   = season_starts(k+1);
        mask = inlet_err_tt.time >= s_start & inlet_err_tt.time < s_end;
        if ~any(mask), continue; end
        err_k = inlet_err_tt.err(mask);
        err_k = err_k(isfinite(err_k));
        if isempty(err_k), continue; end

        mu_k    = mean(err_k);
        sigma_k = std(err_k);
        season_label = sprintf('%d-%02d', year(s_start), mod(year(s_end), 100));
        fprintf('%-25s | %8d | %+9.3f | %8.3f\n', ...
            season_label, numel(err_k), mu_k, sigma_k);

        season_summary = [season_summary; ...
            table({season_label}, s_start, s_end, numel(err_k), mu_k, sigma_k, ...
            'VariableNames', {'season','start','end_','N','mean_err','std_err'})]; %#ok<AGROW>

        fig = figure('Visible', 'off', 'Position', [100 100 800 500]);
        histogram(err_k, 60, 'EdgeColor', 'none');
        grid on;
        xlabel('T_{inlet,fitted} - T_{inlet,ref} [\circC]');
        ylabel('Count');
        title(sprintf('no\\_sensor inlet-temp error, season %s    (N=%d, \\mu=%+.3f, \\sigma=%.3f)', ...
            season_label, numel(err_k), mu_k, sigma_k));
        xline(mu_k, 'r-', 'LineWidth', 1.5);
        xline(mu_k - sigma_k, 'r--');
        xline(mu_k + sigma_k, 'r--');
        save_figure(fig, fullfile(output_folder, sprintf('hist_inlet_err_%s', season_label)));
        close(fig);
    end
    writetable(season_summary, fullfile(output_folder, 'inlet_err_season_summary.csv'));
    fprintf('Season summary + histograms saved to %s\n', output_folder);
end

fprintf('Test complete.\n');
