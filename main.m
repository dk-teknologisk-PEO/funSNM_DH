% main.m
% Main orchestration script for district heating UKF analysis.
% Processes all CSACs simultaneously at each timestep.

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

networks = config.project.datasets.datasets;
output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), '/ukf');
w = waitbar(0.0, "Starting analysis");

% Collect all KPI summaries across all networks
all_kpi_summaries = table();

%% ================================================================
%% MAIN PROCESSING LOOP
%% ================================================================
for network_id = networks
    [meter_data, network_data, topology] = importData(config, network_id);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];
    num_csacs = numel(csac_ids);

    %% ============================================================
    %% INITIALIZE ALL CSACS
    %% ============================================================
    fprintf('\n=== Initializing network %d: %d CSACs ===\n', network_id, num_csacs);
    [all_cs, all_true_traj] = initialize_all_csacs(topology, meter_data, timestamps, ...
        Q_base, R_base, P_base, innovation_gate_initial, config, network_id, drift_config);

    %% ============================================================
    %% TIMESTEP LOOP (all CSACs processed per timestep)
    %% ============================================================
    fprintf('Starting timestep loop: %d timesteps, %d CSACs\n', length(timestamps), num_csacs);

    for t = 1:length(timestamps)
        if mod(t, 50) == 0
            waitbar(t / length(timestamps), w, ...
                sprintf("Timestep %d/%d, network %d (%d CSACs)", ...
                t, length(timestamps), network_id, num_csacs));
        end

        time = timestamps(t);

        % Look up weather data once per timestep (shared across all CSACs)
        current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
        if isempty(current_T_soil_C)
            continue
        end

        % Process each CSAC at this timestep
        for c = 1:num_csacs
            csac_id = csac_ids(c);
            cs = all_cs{c};

            % Extract this CSAC's data for this timestep
            current_data = cs.meter_data(cs.meter_data.timestamp == time, :);
            current_data = sortrows(current_data, 'house_id');

            if isempty(current_data)
                continue
            end

            % Process one timestep for this CSAC
            cs = process_csac_timestep(cs, t, time, current_data, current_T_soil_C, ...
                daily_T_air_max_table, P_base, config, csac_id);

            % Store updated state back
            all_cs{c} = cs;
        end
    end

    %% ============================================================
    %% POST-PROCESSING (all CSACs)
    %% ============================================================
    fprintf('\n=== Post-processing network %d ===\n', network_id);
    network_kpis = post_process_all_csacs(all_cs, all_true_traj, csac_ids, ...
        network_id, output_folder_ukf, kpi_config);
    all_kpi_summaries = [all_kpi_summaries; network_kpis]; %#ok<AGROW>
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