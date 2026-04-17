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

% Load weather data and build lookup table
[T_soil_C, T_air_C] = soilTemp(config);
daily_T_air_max_table = build_daily_T_air_max_table(T_air_C);

networks = config.project.datasets.datasets;
output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), '/ukf');
w = waitbar(0.0, "Starting analysis");

%% ================================================================
%% MAIN PROCESSING LOOP
%% ================================================================
for network = networks
    [meter_data, network_data, topology] = importData(config, network);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];

    for csac = csac_ids

        %% Initialize CSAC state
        houses_csac = ([topology.houses.cul_de_sac_id] == csac);
        topology_csac = topology.houses(houses_csac);
        cs = initialize_csac_state(topology, topology_csac, houses_csac, ...
            meter_data, length(timestamps), Q_base, R_base, P_base, innovation_gate_initial);

        %% Timestep loop
        waitbar(csac / length(csac_ids), w, sprintf("CSAC %d/%d, network %s", ...
            csac, length(csac_ids)-1, string(network)));

        for t = 1:length(timestamps)
            if mod(t, 10) == 0
                waitbar(t / (length(timestamps) * length(csac_ids)) + csac / length(csac_ids), w, ...
                    sprintf("Timestep %d/%d, CSAC %d/%d, network %s", ...
                    t, length(timestamps), csac, length(csac_ids)-1, string(network)));
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

        %% Post-processing: statistics and output
        print_csac_summary(csac, cs.gate_accept_count, cs.gate_reject_count, cs.season_state);

        plot_diagnostics(cs.logger, cs.ground_truth, csac, network, output_folder_ukf);
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

try close(w); catch; end
disp('Analysis complete. All diagnostic plots have been saved.')