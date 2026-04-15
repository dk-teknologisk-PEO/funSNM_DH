% main.m — Main orchestration script
clear all; close all; clc

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', ...
        'src/diagnostics', 'src/particle_filter', 'src/gates', ...
        'src/processing', 'config')

%% 1. Load configuration and weather data
[config, params] = load_config("config.json");
[T_soil_C, T_air_C] = soilTemp(config);
daily_T_air_max_table = precompute_daily_T_air_max(T_air_C);

%% 2. Set up output folders
output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), 'ukf');
output_folder_pf  = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), 'pf');

networks = config.project.datasets.datasets;
w = waitbar(0.0, "Starting analysis");

%% 3. Process each network
for network = networks
    [meter_data, network_data, topology] = importData(config, network);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];

    %% 4. Process each cul-de-sac
    for csac = csac_ids
        % Initialize all states, loggers, season tracking
        csac_state = initialize_csac(topology, meter_data, csac, params, timestamps);

        % Process each timestep
        for t = 1:length(timestamps)
            if mod(t, 10) == 0
                waitbar(t / (length(timestamps) * length(csac_ids)) + csac / length(csac_ids), w, ...
                    sprintf("Timestep %d/%d, CSAC %d, Network %d", t, length(timestamps), csac, network));
            end

            csac_state = process_timestep(csac_state, t, timestamps(t), ...
                T_soil_C, T_air_C, daily_T_air_max_table, params);
        end

        % Save results and print statistics
        print_csac_statistics(csac_state, csac);
        save_all_csac_results(csac_state, csac, network, output_folder_ukf, output_folder_pf);
    end
end

close(w);
disp('Analysis complete. All diagnostic plots have been saved.')