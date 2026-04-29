% main.m
% Main orchestration script for district heating UKF analysis.

% === CLEANUP FROM PREVIOUS RUNS ===
close all force
fclose all;
if exist('w', 'var') && isvalid(w)
    close(w);
end
diary off;
clear all
close all
clc
java.lang.System.gc();

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', ...
    'src/diagnostics', 'src/gates', 'src/CSACs', 'config')

% Read config
config = jsondecode(fileread("config.json"));

% UKF configuration
R_base = config.project.initialization.ukf.measurement_noise^2;
Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, ...
               (config.project.initialization.ukf.process_noise_U)^2]);
P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, ...
               (config.project.initialization.ukf.state_uncertainty_U)^2]);
innovation_gate_initial = config.project.initialization.ukf.state_uncertainty_offset ...
                        * config.project.initialization.innovation_gate_N_sigma;

fprintf('Hard innovation gate set to %.2f °C\n', config.project.initialization.max_innovation_C);

% KPI configuration
kpi_config.offset_tolerance = 0.3;
kpi_config.U_tolerance = 0.02;
kpi_config.convergence_P_offset = 0.5;
kpi_config.convergence_P_U = 0.05;
kpi_config.convergence_hold_days = 14;

% Drift configuration
drift_config.type = 'none';
drift_config.house_index = 1;
drift_config.offset_drift_per_year = 0;
drift_config.step_time = NaT;
drift_config.offset_step = 0;

% Load weather data and build lookup table
[T_soil_C, T_air_C] = soilTemp(config);
daily_T_air_max_table = build_daily_T_air_max_table(T_air_C);

networks = config.project.datasets.datasets(:)';
output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), '/ukf');
w = waitbar(0.0, "Starting analysis");

all_kpi_summaries = table();

%% ================================================================
%% MAIN PROCESSING LOOP
%% ================================================================
for network_id = networks
    [meter_data, network_data, topology] = importData(config, network_id);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];
    num_csacs = numel(csac_ids);

    %% Initialize all CSACs
    fprintf('\n=== Initializing network %d: %d CSACs ===\n', network_id, num_csacs);
    [all_cs, all_true_traj] = initialize_all_csacs(topology, meter_data, timestamps, ...
        Q_base, R_base, P_base, innovation_gate_initial, config, network_id, drift_config);

    %% Initialize shared U_csac
    shared_U_csac = all_cs{1}.U_csac;
    U_csac_true = all_cs{1}.U_csac_true;
    for c = 1:num_csacs
        all_cs{c}.U_csac = shared_U_csac;
    end
    U_csac_update_counter = 0;
    U_csac_history = shared_U_csac;
    U_csac_cfg = config.project.csac_U_estimation;

    %% Initialize main pipe U-value
    sim_cfg = config.project.simulation;
    U_main_true = topology.pipe_parameters.main_pipe.insulation_W_m_K;
    rng(network_id * 10000 + 3e6, 'twister');
    if isfield(sim_cfg, 'U_main_init_mean')
        shared_U_main = sim_cfg.U_main_init_mean + sim_cfg.U_main_init_std * randn();
    else
        shared_U_main = 0.20 + 0.05 * randn();
    end
    shared_U_main = max(0.05, min(0.50, shared_U_main));
    rng('shuffle');
    U_main_update_counter = 0;
    U_main_history = shared_U_main;
    U_main_cfg = config.project.main_pipe_U_estimation;
    main_coupling_cfg = config.project.main_pipe_coupling;
    main_coupling_counter = 0;

    fprintf('  Shared U_csac: init=%.4f, true=%.4f W/m/K\n', shared_U_csac, U_csac_true);
    fprintf('  Shared U_main: init=%.4f, true=%.4f W/m/K\n', shared_U_main, U_main_true);

    %% Timestep loop
    fprintf('Starting timestep loop: %d timesteps, %d CSACs\n', length(timestamps), num_csacs);

    for t = 1:length(timestamps)
        if mod(t, 50) == 0
            waitbar(t / length(timestamps), w, ...
                sprintf("t=%d/%d, net=%d, Uc=%.4f, Um=%.4f", ...
                t, length(timestamps), network_id, shared_U_csac, shared_U_main));
        end

        time = timestamps(t);
        current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
        if isempty(current_T_soil_C)
            continue
        end

%% ============================================================
        %% PROCESS ALL CSACs (single pass)
        %% ============================================================
        any_csac_active = false;
        for c = 1:num_csacs
            cs = all_cs{c};
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

        %% ============================================================
        %% MAIN PIPE COUPLING: Fit and blend for next timestep
        %% ============================================================
        main_coupling_counter = main_coupling_counter + 1;

        if main_coupling_cfg.enabled && main_coupling_counter >= main_coupling_cfg.warmup_timesteps
            [T_csac_inlet_C, main_diag] = estimate_main_pipe_temp(...
                all_cs, csac_ids, topology, current_T_soil_C, shared_U_main, config);

            if main_diag.valid
                % Compute blend alpha based on U_main stability
                % Alpha grows from 0 to max as U_main stabilizes
                if numel(U_main_history) >= 3
                    recent_changes = abs(diff(U_main_history(max(1,end-9):end)));
                    avg_change = mean(recent_changes);
                    % If U_main is changing by less than 0.001 per step, it's stable
                    stability = max(0, 1 - avg_change / 0.002);
                    alpha_blend = min(main_coupling_cfg.max_blend_alpha, ...
                        main_coupling_cfg.max_blend_alpha * stability);
                else
                    alpha_blend = 0;
                end

                % Provide blended T_inlet to each CSAC for next timestep
                for c = 1:num_csacs
                    if isfinite(T_csac_inlet_C(c)) && isfield(all_cs{c}, 'T_inlet_fitted') && ...
                       isfinite(all_cs{c}.T_inlet_fitted)

                        T_blended = (1 - alpha_blend) * all_cs{c}.T_inlet_fitted + ...
                                    alpha_blend * T_csac_inlet_C(c);
                        all_cs{c}.T_inlet_from_main = T_blended;
                    else
                        all_cs{c}.T_inlet_from_main = NaN;
                    end
                    all_cs{c}.main_pipe_blend_alpha = alpha_blend;
                end

                if mod(main_coupling_counter, 100) == 0
                    fprintf('  Main pipe blend: alpha=%.3f, U_main=%.4f\n', ...
                        alpha_blend, shared_U_main);
                end
            end
        end
        
        %% ============================================================
        %% U_csac ESTIMATION (periodically)
        %% ============================================================
        if U_csac_cfg.enabled
            U_csac_update_counter = U_csac_update_counter + 1;
            if U_csac_update_counter >= U_csac_cfg.warmup_timesteps && ...
               mod(U_csac_update_counter, U_csac_cfg.update_interval_timesteps) == 0

                [U_csac_new, U_csac_diag] = update_shared_U_csac(all_cs, shared_U_csac, config);
                if U_csac_diag.adjusted
                    fprintf('  U_csac: %.4f -> %.4f (grad=%.4f, corr=%.2f, houses=%d)\n', ...
                        shared_U_csac, U_csac_new, U_csac_diag.avg_gradient_offset, ...
                        U_csac_diag.avg_corr, U_csac_diag.total_valid_houses);
                end
                shared_U_csac = U_csac_new;
                for c = 1:num_csacs
                    all_cs{c}.U_csac = shared_U_csac;
                end
                U_csac_history(end+1) = shared_U_csac;
            end
        end

        %% ============================================================
        %% U_main ESTIMATION (periodically)
        %% ============================================================
        if U_main_cfg.enabled
            U_main_update_counter = U_main_update_counter + 1;
            if U_main_update_counter >= U_main_cfg.warmup_timesteps && ...
               mod(U_main_update_counter, U_main_cfg.update_interval_timesteps) == 0

                [U_main_new, U_main_diag] = estimate_main_pipe_U(...
                    all_cs, csac_ids, topology, shared_U_main, config);
                if U_main_diag.adjusted
                    fprintf('  U_main: %.4f -> %.4f (grad=%.4f, corr=%.2f)\n', ...
                        shared_U_main, U_main_new, U_main_diag.total_gradient, ...
                        U_main_diag.correlation);
                end
                shared_U_main = U_main_new;
                U_main_history(end+1) = shared_U_main;
            end
        end
    end

    %% ============================================================
    %% CONVERGENCE PLOTS
    %% ============================================================
    if ~exist(output_folder_ukf, 'dir'), mkdir(output_folder_ukf); end

    % U_csac convergence
    fig_Uc = figure('Visible', 'off', 'Position', [100, 100, 800, 400]);
    plot(1:numel(U_csac_history), U_csac_history, 'b-', 'LineWidth', 2, 'DisplayName', 'Estimated');
    yline(U_csac_true, 'r--', 'LineWidth', 2, 'DisplayName', 'True');
    xlabel('Update step'); ylabel('U_{csac} [W/m/K]');
    title(sprintf('Network %d — U_{csac} (init=%.4f, final=%.4f, true=%.4f)', ...
        network_id, U_csac_history(1), U_csac_history(end), U_csac_true));
    legend('Location', 'best'); grid on;
    save_figure(fig_Uc, fullfile(output_folder_ukf, sprintf('U_csac_convergence_network_%d', network_id)));
    close(fig_Uc);

    % U_main convergence
    fig_Um = figure('Visible', 'off', 'Position', [100, 100, 800, 400]);
    plot(1:numel(U_main_history), U_main_history, 'b-', 'LineWidth', 2, 'DisplayName', 'Estimated');
    yline(U_main_true, 'r--', 'LineWidth', 2, 'DisplayName', 'True');
    xlabel('Update step'); ylabel('U_{main} [W/m/K]');
    title(sprintf('Network %d — U_{main} (init=%.4f, final=%.4f, true=%.4f)', ...
        network_id, U_main_history(1), U_main_history(end), U_main_true));
    legend('Location', 'best'); grid on;
    save_figure(fig_Um, fullfile(output_folder_ukf, sprintf('U_main_convergence_network_%d', network_id)));
    close(fig_Um);

    %% ============================================================
    %% POST-PROCESSING
    %% ============================================================
    fprintf('\n=== Post-processing network %d ===\n', network_id);
    fprintf('  Final U_csac: %.4f (true=%.4f, error=%.4f)\n', ...
        shared_U_csac, U_csac_true, shared_U_csac - U_csac_true);
    fprintf('  Final U_main: %.4f (true=%.4f, error=%.4f)\n', ...
        shared_U_main, U_main_true, shared_U_main - U_main_true);

    plot_csac_U_diagnostics(all_cs, csac_ids, shared_U_csac, network_id, output_folder_ukf);

    network_kpis = post_process_all_csacs(all_cs, all_true_traj, csac_ids, ...
        network_id, output_folder_ukf, kpi_config);
    all_kpi_summaries = [all_kpi_summaries; network_kpis]; %#ok<AGROW>
end

%% Save aggregate KPI summary
if ~isempty(all_kpi_summaries)
    if ~exist(output_folder_ukf, 'dir'), mkdir(output_folder_ukf); end
    writetable(all_kpi_summaries, fullfile(output_folder_ukf, 'kpi_summary_all.csv'));
    fprintf('\nAggregate KPI summary saved: %d houses across %d CSACs\n', ...
        height(all_kpi_summaries), numel(unique(all_kpi_summaries.csac_id)));
end

try close(w); catch; end
disp('Analysis complete. All diagnostic plots have been saved.')