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
kpi_config.offset_tolerance = 0.3;       % °C
kpi_config.U_tolerance = 0.02;           % W/m/K
kpi_config.convergence_P_offset = 0.5;   % °C std threshold
kpi_config.convergence_P_U = 0.05;       % W/m/K std threshold
kpi_config.convergence_hold_days = 14;   % active days

% Drift configuration — change these to test different scenarios
drift_config.type = 'none';              % 'none', 'linear', or 'step'
drift_config.offset_drift_per_year = 0;  % °C/year (for 'linear')
drift_config.step_time = NaT;            % datetime (for 'step')
drift_config.offset_step = 0;            % °C (for 'step')

% Load weather data and build lookup table
[T_soil_C, T_air_C] = soilTemp(config);
daily_T_air_max_table = build_daily_T_air_max_table(T_air_C);

networks = config.project.datasets.datasets;
output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), '/ukf');
w = waitbar(0.0, "Starting analysis");

% Collect all KPI summaries across the run
all_kpi_summaries = table();

%% ================================================================
%% MAIN PROCESSING LOOP
%% ================================================================
for network_id = networks
    [meter_data, network_data, topology] = importData(config, network_id);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];

    for csac = csac_ids

        %% Initialize CSAC state (includes base offset generation and application)
        houses_csac = ([topology.houses.cul_de_sac_id] == csac);
        topology_csac = topology.houses(houses_csac);
        cs = initialize_csac_state(topology, topology_csac, houses_csac, ...
            meter_data, length(timestamps), Q_base, R_base, P_base, ...
            innovation_gate_initial, config, network_id);

        %% Apply offset drift to meter data (additional to base offset)
        cs.meter_data = apply_offset_drift_to_data(...
            cs.meter_data, cs.house_ids, drift_config, timestamps);

        %% Generate true trajectories for KPI evaluation
        true_trajectories = generate_true_trajectories(...
            cs.ground_truth, timestamps, drift_config);

        %% Timestep loop
        waitbar(csac / length(csac_ids), w, sprintf("CSAC %d/%d, network %d", ...
            csac, length(csac_ids)-1, network_id));

        for t = 1:length(timestamps)
            if mod(t, 10) == 0
                waitbar(t / (length(timestamps) * length(csac_ids)) + csac / length(csac_ids), w, ...
                    sprintf("Timestep %d/%d, CSAC %d/%d, network %d", ...
                    t, length(timestamps), csac, length(csac_ids)-1, network_id));
            end

            time = timestamps(t);
            current_data = cs.meter_data(cs.meter_data.timestamp == time, :);
            current_data = sortrows(current_data, 'house_id');
            current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;

            if isempty(current_data) || isempty(current_T_soil_C)
                continue
            end

            cs = process_csac_timestep(cs, t, time, current_data, current_T_soil_C, ...
                daily_T_air_max_table, P_base, config, csac);
        end

        %% Post-processing: statistics, KPIs, and output
        print_csac_summary(csac, cs.gate_accept_count, cs.gate_reject_count, cs.season_state);

        kpi_summary = compute_and_save_network_kpis(cs, csac, output_folder_ukf, ...
            kpi_config, true_trajectories);
        all_kpi_summaries = [all_kpi_summaries; kpi_summary]; %#ok<AGROW>

        plot_diagnostics(cs.logger, cs.ground_truth, csac, network_id, output_folder_ukf);
        save_logger_to_csv(cs.logger, output_folder_ukf, strcat('ukf_csac_', string(csac)));
        save_diagnostic_summary(cs.logger, cs.ground_truth, csac, output_folder_ukf, 'ukf');
        save_diagnostic_summary_detailed(cs.logger, cs.ground_truth, csac, output_folder_ukf, 'ukf');
        save_daily_diagnostics(cs.logger, cs.ground_truth, csac, output_folder_ukf, 'ukf');

        figs = findall(0, 'Type', 'figure');
        for fig = figs'
            if ~strcmp(get(fig, 'Tag'), 'TMWWaitbar')
                close(fig);
            end
        end
    end
end

%% Save aggregate KPI summary
if ~isempty(all_kpi_summaries)
    if ~exist(output_folder_ukf, 'dir')
        mkdir(output_folder_ukf);
    end
    writetable(all_kpi_summaries, fullfile(output_folder_ukf, 'kpi_summary_all.csv'));
    fprintf('\nAggregate KPI summary saved: %d houses across %d CSACs\n', ...
        height(all_kpi_summaries), numel(unique(all_kpi_summaries.csac_id)));
end

try close(w); catch; end
disp('Analysis complete. All diagnostic plots have been saved.')