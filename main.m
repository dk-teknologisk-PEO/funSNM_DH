% main.m
% Clear workspace, add paths, set parameters

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

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'src/gates', 'src/CSACs/', 'config')

% read config-file
config = jsondecode(fileread("config.json"));

% --- DEBUG SETTINGS (from config) ---
debug_disable_stability_gate = config.project.debug.disable_stability_gate;
debug_print_update_result = config.project.debug.print_update_result;

% --- HARD INNOVATION GATE (from config) ---
hard_innovation_gate_C = config.project.initialization.max_innovation_C;
max_nis = config.project.initialization.max_nis;
fprintf('Hard innovation gate set to %.2f °C\n', hard_innovation_gate_C);

% ukf-configuration
R_base = config.project.initialization.ukf.measurement_noise^2;
Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, (config.project.initialization.ukf.process_noise_U)^2]);
P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, (config.project.initialization.ukf.state_uncertainty_U)^2]);
innovation_gate_initial = config.project.initialization.ukf.state_uncertainty_offset * config.project.initialization.innovation_gate_N_sigma;

% gating configuration
absolute_flow_floor_kg_h = config.project.cutoff.flow_cutoff;
delta_T_gate_threshold = config.project.cutoff.delta_T_gate_threshold;
alpha_min = config.project.cutoff.alpha_min;
max_delta_T_change = config.project.initialization.max_delta_T_change_rate;

networks = config.project.datasets.datasets;

% load network topology and measurement data
[T_soil_C, T_air_C] = soilTemp(config);

%% ============================================================
%% PRE-COMPUTE DAILY T_AIR_MAX LOOKUP TABLE
%% ============================================================
fprintf('Pre-computing daily T_air_max lookup table...\n');

T_air_C.date = dateshift(T_air_C.time, 'start', 'day');
[air_date_groups, air_unique_dates] = findgroups(T_air_C.date);

daily_T_air_max_values = splitapply(@max, T_air_C.values, air_date_groups);
daily_T_air_max_table = table(air_unique_dates, daily_T_air_max_values, ...
    'VariableNames', {'date', 'T_air_max'});
daily_T_air_max_table = sortrows(daily_T_air_max_table, 'date');

fprintf('  Built daily T_air_max table: %d days (%.1f to %.1f °C)\n', ...
    height(daily_T_air_max_table), ...
    min(daily_T_air_max_table.T_air_max), ...
    max(daily_T_air_max_table.T_air_max));

output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), '/ukf');
w = waitbar(0.0, "Starting analysis");

% run through all the networks from the config-file
for network = networks
    % load relevant data for the network
    [meter_data, network_data, topology] = importData(config, network);
    timestamps = unique(meter_data.timestamp);
    
    % find the ids of the csacs from the topology
    csac_ids = [topology.cul_de_sacs.id];

    % loop through all cul-de-sacs in the network
    for csac = csac_ids
        
        %% ============================================================
        %% CSAC INITIALIZATION
        %% ============================================================
        houses_csac = ([topology.houses.cul_de_sac_id] == csac);
        topology_csac = topology.houses(houses_csac);

        cs = initialize_csac_state(topology, topology_csac, houses_csac, ...
            meter_data, length(timestamps), Q_base, R_base, P_base, innovation_gate_initial);

        % run through the timesteps
        waitbar((csac)/length(csac_ids), w, strcat("Running network. At timestep: 0/", ...
            string(length(timestamps)), ", csac: ", string(csac), "/", ...
            string(length(csac_ids)-1), ", network: ", string(network)))

        for t = 1:length(timestamps)

            if mod(t, 10) == 0
                waitbar(t/(length(timestamps)*length(csac_ids)) + (csac)/length(csac_ids), w, ...
                    strcat("Running network. At timestep: ", string(t), "/", ...
                    string(length(timestamps)), ", csac: ", string(csac), "/", ...
                    string(length(csac_ids)-1), ", network: ", string(network)))
            end
            time = timestamps(t);

            current_data = cs.meter_data(cs.meter_data.timestamp == time, :);
            current_data = sortrows(current_data, 'house_id');

            current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;

            if isempty(current_data) || isempty(current_T_soil_C)
                continue
            end

            %% ============================================================
            %% HEATING SEASON GATE
            %% ============================================================
            current_date = dateshift(time, 'start', 'day');
            [cs.season_state, actions] = manage_heating_season(...
                cs.season_state, current_date, t, daily_T_air_max_table, config, csac);

            if actions.do_restore_snapshot || actions.do_save_snapshot
                [cs.ukf_states, cs.ukf_states_snapshot] = apply_season_actions(...
                    actions, cs.ukf_states, cs.ukf_states_snapshot, P_base, config);
            end

            if actions.season_is_active
                cs.season_state.last_active_date = current_date;
            end

            if actions.skip_timestep
                if t > 1
                    for i = 1:cs.num_houses
                        cs.logger.state_estimates(:, i, t) = cs.logger.state_estimates(:, i, t-1);
                        cs.logger.covariance_posterior(:, i, t) = cs.logger.covariance_posterior(:, i, t-1);
                    end
                end
                cs.logger.timestamps(t) = time;
                continue;
            end

            %% ============================================================
            %% ACTIVE SEASON: Main pipe temperature estimation
            %% ============================================================
            num_active_houses = sum(current_data.flow_kg_h >= absolute_flow_floor_kg_h);
            is_csac_active = (num_active_houses >= config.project.initialization.min_active_houses);

            [T_junction_ukf_C, ukf_master_offset] = calculate_main_pipe_temp(...
                current_data, current_T_soil_C, cs.U_csac, absolute_flow_floor_kg_h, cs.ukf_states);
            
            if all(isnan(T_junction_ukf_C))
                continue
            end

            current_data.T_main_ukf_C = T_junction_ukf_C;
            cs.logger.timestamps(t) = time;

            %% ============================================================
            %% PER-HOUSE UKF UPDATES
            %% ============================================================
            for i = 1:cs.num_houses
                log_ukf = false;
                house_id = cs.house_ids(i);
                house_data = current_data(current_data.house_id == house_id, :);
                if isempty(house_data)
                    continue
                end

                % --- Flow and temperature validity gates ---
                [is_flow_sufficient, ~, ~, ~] = flow_validity_gate(...
                    house_data.flow_kg_h, cs.ukf_states{i}.x(2), house_data.length_service_m, alpha_min);
                is_flow_sufficient = is_flow_sufficient && (house_data.flow_kg_h > absolute_flow_floor_kg_h);
                delta_T_sufficient = ((house_data.T_main_ukf_C - current_T_soil_C) >= delta_T_gate_threshold);
                can_update = is_csac_active && is_flow_sufficient && delta_T_sufficient;

                % --- Gated UKF update ---
                if can_update
                    gate_params = struct( ...
                        'last_valid_T_main_C',    cs.last_valid_T_main_C(i), ...
                        'max_delta_T_change',     max_delta_T_change, ...
                        'hard_innovation_gate_C', hard_innovation_gate_C, ...
                        'max_nis',                max_nis, ...
                        'disable_stability_gate', debug_disable_stability_gate, ...
                        'debug_print',            debug_print_update_result, ...
                        'house_id',               house_id, ...
                        'time',                   time);

                    [cs.ukf_states{i}, update_result] = update_house_ukf_gated(...
                        cs.ukf_states{i}, house_data, current_T_soil_C, config, gate_params);

                    if update_result.accepted
                        cs.gate_accept_count = cs.gate_accept_count + 1;
                        cs.ukf_offsets(i) = cs.ukf_states{i}.x(1);
                        log_ukf = true;
                        cs.last_valid_T_main_C(i) = house_data.T_main_ukf_C;
                        cs.last_update_timestamp(i) = time;
                        if isfinite(update_result.innovation_gate)
                            cs.ukf_innovation_gate(i) = update_result.innovation_gate;
                        end
                    elseif update_result.rejected
                        cs.gate_reject_count = cs.gate_reject_count + 1;
                    end
                end

                % --- Logging ---
                if log_ukf
                    cs.logger = update_logger(cs.logger, t, i, time, cs.ukf_states{i}, update_result.diagnostics);
                elseif t > 1
                    cs.logger.state_estimates(:, i, t) = cs.logger.state_estimates(:, i, t-1);
                    cs.logger.covariance_posterior(:, i, t) = cs.logger.covariance_posterior(:, i, t-1);
                end
            end

            %% ============================================================
            %% SNAPSHOT UPDATE
            %% ============================================================
            [cs.ukf_states_snapshot, cs.season_state.snapshot_timestep] = ...
                update_snapshot(cs.season_state, current_date, t, ...
                cs.ukf_states, cs.ukf_states_snapshot, config);

            %% ============================================================
            %% MASTER OFFSET APPLICATION
            %% ============================================================
            cs.ukf_states = apply_master_offset(cs.ukf_states, ukf_master_offset, ...
                num_active_houses, config);
        end

        %% ============================================================
        %% POST-PROCESSING: Statistics and Output
        %% ============================================================
        print_csac_summary(csac, cs.gate_accept_count, cs.gate_reject_count, cs.season_state);

        plot_diagnostics(cs.logger, cs.ground_truth, csac, network, output_folder_ukf);
        save_logger_to_csv(cs.logger, output_folder_ukf, strcat('ukf_csac_', string(csac)));
        save_diagnostic_summary(cs.logger, cs.ground_truth, csac, output_folder_ukf, 'ukf');
        save_diagnostic_summary_detailed(cs.logger, cs.ground_truth, csac, output_folder_ukf, 'ukf');
        save_daily_diagnostics(cs.logger, cs.ground_truth, csac, output_folder_ukf, 'ukf');

        % Close diagnostic figures but keep waitbar
        figs = findall(0, 'Type', 'figure');
        for fig = figs'
            if ~strcmp(get(fig, 'Tag'), 'TMWWaitbar')
                close(fig);
            end
        end
    end
end

try
    close(w);
catch
end
disp('Analysis complete. All diagnostic plots have been saved.')