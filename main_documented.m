% =========================================================================
% Main script for estimating district heating service pipe insulation (U-value)
% and utility meter temperature offsets using Unscented Kalman Filters (UKF)
% and Particle Filters (PF).
%
% The script processes data for specified networks and cul-de-sacs,
% applying a sophisticated gating mechanism to ensure data quality before
% updating the filter states.
%
% Features:
% - Dual filtering approach (UKF and PF) for comparison.
% - Dynamic N-sigma gating for innovations.
% - Transient detection based on temperature change rate.
% - State hibernation and uncertainty re-inflation after long data gaps.
% - Efficient two-tier temperature calculation to improve performance.
% - Automated plotting and CSV export of results.
% =========================================================================

%% --- 1. Initialization and Configuration ---
clear all;
close all;
clc;

% Add function paths
addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'config')

% Load configuration from JSON file
config = jsondecode(fileread("config.json"));

% --- Filter Parameters ---
R_base = config.project.initialization.ukf.measurement_noise^2;
Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, (config.project.initialization.ukf.process_noise_U)^2]);
P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, (config.project.initialization.ukf.state_uncertainty_U)^2]);
num_particles = config.project.initialization.pf.num_particles;

% --- Gating and Logic Parameters ---
flow_threshold = config.project.cutoff.flow_cutoff;
delta_T_gate_threshold = config.project.cutoff.delta_T_gate_threshold;
min_active_houses = config.project.initialization.min_active_houses;
gate_N_sigma = config.project.initialization.innovation_gate_N_sigma;
max_delta_T_change_rate = config.project.initialization.max_delta_T_change_rate;
air_temp_cutoff = config.project.initialization.max_air_temperature;
hibernation_threshold_h = config.project.initialization.hibernation_reset_threshold_hours;
error_meaning = config.project.initialization.modify_mean_offset;

networks = config.project.datasets.datasets;

% --- Setup Output and Data ---
[T_soil_C, T_air_C] = soilTemp(config);
output_folder = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'));
output_folder_ukf = fullfile(output_folder, 'ukf');
output_folder_pf = fullfile(output_folder, 'pf');
w = waitbar(0.0, "Starting analysis");

%% --- 2. Main Processing Loop ---
for network = networks
    % --- Network Data Loading ---
    [meter_data, network_data, topology] = importData(config, network);
    timestamps = unique(meter_data.timestamp);
    csac_ids = [topology.cul_de_sacs.id];

    % --- Cul-de-Sac Loop ---
    for csac = csac_ids
        
        % --- CSAC-Specific Setup ---
        U_csac = topology.pipe_parameters.csac_pipe.insulation_W_m_K;
        houses_csac_mask = ([topology.houses.cul_de_sac_id] == csac);
        topology_csac = topology.houses(houses_csac_mask);
        num_houses_csac = sum(houses_csac_mask);
        houses_csac_ids = [topology_csac.id];

        % Prepare ground truth and geometry tables for this CSAC
        ground_truth_csac = table([topology_csac.id]', zeros(num_houses_csac, 1), [topology_csac.service_pipe_insulation_W_m_K]', ...
            'VariableNames', {'house_id', 'true_offset', 'true_U'});
        csac_table = table([topology_csac.id]', [topology_csac.service_pipe_len_m]', [topology_csac.dist_on_cul_de_sac_m]', ...
            'VariableNames', {'house_id', 'length_service_m', 'x_pos_m'});
        meter_data_csac = join(meter_data(ismember(meter_data.house_id, houses_csac_ids), :), csac_table, "Keys", "house_id");

        % --- Filter and State Initialization for this CSAC ---
        ukf_states = cell(1, num_houses_csac);
        pf_particles = cell(1, num_houses_csac);
        pf_states = cell(1, num_houses_csac);
        for i = 1:num_houses_csac
            x_init = [randn()*0.7; 0.12 + randn()*0.03]; % Randomized initial guess
            ukf_states{i} = initialize_kalman_filter(x_init, Q_base, R_base, P_base);
            pf_particles{i} = initialize_pf_state(num_particles, x_init, P_base);
            pf_states{i}.x = x_init;
            pf_states{i}.P = P_base;
        end
        
        % --- Logger and Tracker Initialization ---
        logger_ukf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);
        logger_pf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);
        ukf_offsets = nan(num_houses_csac, 1);
        pf_offsets = nan(num_houses_csac, 1);
        
        % State trackers for advanced gating and hibernation
        last_valid_T_main_ukf_C = nan(num_houses_csac, 1);
        last_valid_T_main_pf_C = nan(num_houses_csac, 1);
        last_update_timestamp_ukf = NaT(num_houses_csac, 1);
        last_update_timestamp_pf = NaT(num_houses_csac, 1);
        waitbar_msg = sprintf("Network %s, CSAC %d, Timestep %d/%d", string(network), csac, 0, length(timestamps));
        waitbar(0, w, waitbar_msg);
        % --- Timestep Loop ---
        for t = 1:length(timestamps)
            if mod(t, 10) == 0
                waitbar_msg = sprintf("Network %s, CSAC %d, Timestep %d/%d", string(network), csac, t, length(timestamps));
                waitbar(t/length(timestamps), w, waitbar_msg);
            end
            
            time = timestamps(t);
            current_data = meter_data_csac(meter_data_csac.timestamp == time, :);
            
            % --- Timestep Validity Check ---
            if isempty(current_data)
                continue; % Skip if no data for any house at this time
            end
            current_T_soil_C = T_soil_C.values(T_soil_C.time == time);
            current_T_air_C = T_air_C.values(year(T_air_C.time)==year(time) & month(T_air_C.time)==month(time) & day(T_air_C.time)==day(time));
            if isempty(current_T_soil_C) || isempty(current_T_air_C)
                continue; % Skip if weather data is missing or too warm
            end

            T_junction_ukf_C = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_threshold, ukf_states);
            current_data.T_main_ukf_C = T_junction_ukf_C;

            T_junction_pf_C = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_threshold, pf_states);
            current_data.T_main_pf_C = T_junction_pf_C;
            
            % --- Tier 1: Fast, Approximate Temperature Calculation ---
            % Used for initial gating checks for all houses.
            % u_service_ukf = cellfun(@(s) s.x(2), ukf_states);
            % offsets_ukf   = cellfun(@(s) s.x(1), ukf_states);
            % current_data.T_main_approx_ukf_C = get_main_temp(current_data.T_supply_C + offsets_ukf', current_data.flow_kg_h, u_service_ukf', current_data.length_service_m, current_T_soil_C);
            % 
            % u_service_pf = cellfun(@(s) s.x(2), pf_states);
            % offsets_pf   = cellfun(@(s) s.x(1), pf_states);
            % current_data.T_main_approx_pf_C = get_main_temp(current_data.T_supply_C + offsets_pf', current_data.flow_kg_h, u_service_pf', current_data.length_service_m, current_T_soil_C);

            is_csac_active = (sum(current_data.flow_kg_h >= flow_threshold) > min_active_houses);

            % --- House Loop ---
            for i = 1:num_houses_csac
                log_ukf = false;
                log_pf = false;
                house_data = current_data(current_data.house_id == houses_csac_ids(i), :);
                if isempty(house_data)
                    continue; 
                end

                % ======================= UKF Update Logic =======================
                % 1. Primary Gating (Fast Checks)
                is_flow_sufficient = (house_data.flow_kg_h >= flow_threshold);
                is_delta_T_sufficient = ((house_data.T_main_ukf_C - current_T_soil_C) >= delta_T_gate_threshold);
                if is_csac_active && is_flow_sufficient && is_delta_T_sufficient

                    % 2. Stability Gating (Transient Rejection)
                    if isnan(last_valid_T_main_ukf_C(i)) 
                        delta_T_main_change = 0;
                    else
                        delta_T_main_change = abs(house_data.T_main_ukf_C - last_valid_T_main_ukf_C(i)); 
                    end
                    
                    last_valid_T_main_ukf_C(i) = house_data.T_main_ukf_C; % Update baseline
                    
                    if delta_T_main_change < max_delta_T_change_rate
                        
                        % 3. Innovation Gating (Dynamic Outlier Rejection)
                        % Tier 2: Precise temperature calculation, only run when needed.
                        % T_junction_vector = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_threshold, ukf_states);
                        % house_data.T_main_ukf_C = T_junction_vector(i);
                        
                        [predicted_temp, innovation_var] = predict_measurement_ukf(ukf_states{i}, house_data, current_T_soil_C);
                        innovation = house_data.T_supply_C - predicted_temp;
                        innovation_gate = gate_N_sigma * sqrt(innovation_var);

                        if abs(innovation) < innovation_gate
                            % 4. Hibernation Check (Reset uncertainty after long gap)
                            time_gap_h = hours(time - last_update_timestamp_ukf(i));
                            if ~isnat(last_update_timestamp_ukf(i)) && time_gap_h > hibernation_threshold_h
                                fprintf('UKF Hibernation Reset for House %d at %s (Gap: %.1f hours)\n', houses_csac_ids(i), datestr(time), time_gap_h);
                                ukf_states{i}.P(1,1) = P_base(1,1); % Reset offset uncertainty
                                ukf_states{i}.P(2,2) = ukf_states{i}.P(2,2) * 1.5; % Slightly increase U-value uncertainty
                                last_valid_T_main_ukf_C(i) = NaN; % Reset stability tracker
                            end
                            
                            % 5. Execute Filter Update
                            [ukf_states{i}, diagnostics_ukf] = update_ukf_house(ukf_states{i}, house_data, current_T_soil_C, config);
                            ukf_offsets(i) = ukf_states{i}.x(1);
                            log_ukf = true;
                            last_update_timestamp_ukf(i) = time; % Record timestamp of successful update
                        end
                    end
                end

                % ======================= PF Update Logic ========================
                % (This logic mirrors the UKF structure exactly)
                is_flow_sufficient = (house_data.flow_kg_h >= flow_threshold);
                is_delta_T_sufficient = ((house_data.T_main_pf_C - current_T_soil_C) >= delta_T_gate_threshold);
                if is_csac_active && is_flow_sufficient && is_delta_T_sufficient
                    
                    if isnan(last_valid_T_main_pf_C(i)) 
                        delta_T_main_change = 0;
                    else 
                        delta_T_main_change = abs(house_data.T_main_pf_C - last_valid_T_main_pf_C(i)); 
                    end
                    
                    last_valid_T_main_pf_C(i) = house_data.T_main_pf_C;

                    if delta_T_main_change < max_delta_T_change_rate
                        % T_junction_vector = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_threshold, pf_states);
                        % house_data.T_main_pf_C = T_junction_vector(i);

                        [predicted_temp, innovation_var, prop_particles] = predict_measurement_pf(pf_particles{i}, house_data, current_T_soil_C, R_base, Q_base, config);
                        innovation = house_data.T_supply_C - predicted_temp;
                        innovation_gate = gate_N_sigma * sqrt(innovation_var);

                        if abs(innovation) < innovation_gate
                            time_gap_h = hours(time - last_update_timestamp_pf(i));
                            if ~isnat(last_update_timestamp_pf(i)) && time_gap_h > hibernation_threshold_h
                                fprintf('PF Hibernation Reset for House %d at %s (Gap: %.1f hours)\n', houses_csac_ids(i), datestr(time), time_gap_h);
                                current_mean = pf_states{i}.x;
                                spread_noise = chol(P_base)' * randn(2, num_particles);
                                pf_particles{i} = current_mean + spread_noise;
                                % Re-apply constraints...
                                last_valid_T_main_pf_C(i) = NaN;
                            end
                            
                            [pf_particles{i}, est_pf, cov_pf, diagnostics_pf] = update_pf_house(prop_particles, house_data, current_T_soil_C, R_base, predicted_temp);
                            pf_states{i}.x = est_pf;
                            pf_states{i}.P = cov_pf;
                            pf_offsets(i) = pf_states{i}.x(1);
                            log_pf = true;
                            last_update_timestamp_pf(i) = time;
                        end
                    end
                end

                % --- Update Loggers ---
                if log_ukf
                    logger_ukf = update_logger(logger_ukf, t, i, time, ukf_states{i}, diagnostics_ukf);
                elseif t > 1
                    logger_ukf.state_estimates(:, i, t) = logger_ukf.state_estimates(:, i, t - 1);
                    logger_ukf.covariance_posterior(:, i, t) = logger_ukf.covariance_posterior(:, i, t - 1);
                end
                if log_pf
                    logger_pf = update_logger(logger_pf, t, i, time, pf_states{i}, diagnostics_pf);
                elseif t > 1
                    logger_pf.state_estimates(:, i, t) = logger_pf.state_estimates(:, i, t - 1);
                    logger_pf.covariance_posterior(:, i, t) = logger_pf.covariance_posterior(:, i, t - 1);
                end
            end % End of house loop

            % --- Mean Offset Correction ---
            if error_meaning
                mean_ukf_offsets = mean(ukf_offsets, 'omitmissing');
                mean_pf_offsets = mean(pf_offsets, 'omitmissing');
                for i = 1:num_houses_csac
                    if ~isnan(mean_ukf_offsets) && ~isnan(ukf_states{i}.x(1))
                        ukf_states{i}.x(1) = ukf_states{i}.x(1) - mean_ukf_offsets;
                    end
                    if ~isnan(mean_pf_offsets) && ~isnan(pf_states{i}.x(1))
                        pf_states{i}.x(1) = pf_states{i}.x(1) - mean_pf_offsets;
                        pf_particles{i}(1, :) = pf_particles{i}(1, :) - mean_pf_offsets;
                    end
                end
            end
        end % End of timestep loop

        % --- Save Results for this CSAC ---
        plot_diagnostics(logger_ukf, ground_truth_csac, csac, network, output_folder_ukf);
        save_logger_to_csv(logger_ukf, output_folder_ukf, sprintf('ukf_csac%d', csac));
        plot_diagnostics(logger_pf, ground_truth_csac, csac, network, output_folder_pf);
        save_logger_to_csv(logger_pf, output_folder_pf, sprintf('pf_csac%d', csac));
        
    end % End of CSAC loop
end % End of network loop

close(w);
disp('Analysis complete. All diagnostic plots and CSV files have been saved.')