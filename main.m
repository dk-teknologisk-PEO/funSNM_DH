% Clear workspace, add paths, set parameters
clear all
close all
clc

% --- DEBUG SETTINGS ---
debug_disable_stability_gate = false;
debug_disable_mean_centering = true;
debug_print_update_result = true;


addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'src/gates', 'config')

% read config-file
config = jsondecode(fileread("config.json"));

% ukf-configuration
R_base = config.project.initialization.ukf.measurement_noise^2;
Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, ...
               (config.project.initialization.ukf.process_noise_U)^2]);
P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, ...
               (config.project.initialization.ukf.state_uncertainty_U)^2]);
innovation_gate_initial = config.project.initialization.ukf.state_uncertainty_offset ...
                        * config.project.initialization.innovation_gate_N_sigma;

% pf-configuration
num_particles = config.project.initialization.pf.num_particles;

% gating configuration
absolute_flow_floor_kg_h = config.project.cutoff.flow_cutoff;
delta_T_gate_threshold = config.project.cutoff.delta_T_gate_threshold;
min_active_houses = config.project.initialization.min_active_houses;
gate_N_sigma = config.project.initialization.innovation_gate_N_sigma;
max_innovation = config.project.initialization.max_innovation_C;
alpha_min = config.project.cutoff.alpha_min;

max_delta_T_change = config.project.initialization.max_delta_T_change_rate;

air_temp_cutoff = config.project.initialization.max_air_temperature;

error_meaning = config.project.initialization.modify_mean_offset;
hibernation_reset_threshold_hours = config.project.initialization.hibernation_reset_threshold_hours;

% Master offset damping parameters
master_offset_lr = config.project.initialization.master_offset_learning_rate;
master_offset_threshold = config.project.initialization.master_offset_min_threshold;

networks = config.project.datasets.datasets;

% load network topology and measurement data
[T_soil_C, T_air_C] = soilTemp(config);
output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), 'ukf');
output_folder_pf  = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'), 'pf');
w = waitbar(0.0, "Starting analysis");


% run through all the networks from the config-file
for network = networks
    % load relevant data for the network
    [meter_data, network_data, topology] = importData(config, network);
    timestamps = unique(meter_data.timestamp);
    
    % find the ids of the csacs from the topology
    csac_ids = [topology.cul_de_sacs.id];
    prev_houses = 0;

    % loop through all cul-de-sacs in the network
    for csac = csac_ids
        
        % extract relevant data for the csac
        U_csac = topology.pipe_parameters.csac_pipe.insulation_W_m_K;
        houses_csac = ([topology.houses.cul_de_sac_id] == csac);
        topology_csac = topology.houses(houses_csac);
        num_houses_csac = sum(houses_csac);
        houses_csac_ids = sort([topology.houses(houses_csac).id]);
        true_offset = zeros([size(topology_csac), 1]);

        ground_truth_csac = table([topology_csac.id]', true_offset, ...
            [topology_csac.service_pipe_insulation_W_m_K]', ...
            'VariableNames', {'house_id', 'true_offset', 'true_U'});
        
        csac_table = table([topology_csac.id]', [topology_csac.service_pipe_len_m]', ...
            [topology_csac.dist_on_cul_de_sac_m]', ...
            'VariableNames', {'house_id', 'length_service_m', 'x_pos_m'});

        meter_data_csac = meter_data(ismember(meter_data.house_id, houses_csac_ids), :);
        meter_data_csac = join(meter_data_csac, csac_table, "Keys", "house_id");

        % --- Initialize filters for each house ---
        ukf_states = cell(1, num_houses_csac);
        pf_particles = cell(1, num_houses_csac);
        pf_states = cell(1, num_houses_csac);
        ukf_innovation_gate = nan([1, num_houses_csac]);
        pf_innovation_gate  = nan([1, num_houses_csac]);
        
        last_valid_T_main_ukf_C = nan(num_houses_csac, 1);
        last_valid_T_main_pf_C  = nan(num_houses_csac, 1);
        last_update_timestamp_ukf = NaT(num_houses_csac, 1);
        last_update_timestamp_pf  = NaT(num_houses_csac, 1);

        for i = 1:num_houses_csac
            x_init = [randn()*0.3; 0.12 + randn()*0.03];
            ukf_states{i} = initialize_kalman_filter(x_init, Q_base, R_base, P_base);
            pf_particles{i} = initialize_pf_state(num_particles, x_init, P_base, config);
            pf_states{i}.x = x_init;
            pf_states{i}.P = P_base;
            ukf_innovation_gate(i) = innovation_gate_initial;
            pf_innovation_gate(i)  = innovation_gate_initial;
        end
        
        logger_ukf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);
        logger_pf  = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);

        ukf_offsets = nan([num_houses_csac, 1]);
        pf_offsets  = nan([num_houses_csac, 1]);

        % ================================================================
        % MAIN TIME LOOP
        % ================================================================
        waitbar((csac)/length(csac_ids), w, strcat("Running network. At timestep: 0/", ...
            string(length(timestamps)), ", csac: ", string(csac), "/", ...
            string(length(csac_ids)), ", network: ", string(network)))

        for t = 1:length(timestamps)

            disp(t)
            if mod(t, 10) == 0
                waitbar(t/(length(timestamps)*length(csac_ids)) + (csac)/length(csac_ids), w, ...
                    strcat("Running network. At timestep: ", string(t), "/", ...
                    string(length(timestamps)), ", csac: ", string(csac), "/", ...
                    string(length(csac_ids)), ", network: ", string(network)))
            end
            time = timestamps(t);

            current_data = meter_data_csac(meter_data_csac.timestamp == time, :);
            current_data = sortrows(current_data, 'house_id');

            current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
            current_T_air_C  = T_air_C((year(T_air_C.time) == year(time)) & ...
                (month(T_air_C.time) == month(time)) & ...
                (day(T_air_C.time) == day(time)), :).values;

            if isempty(current_data) || isempty(current_T_soil_C)
                continue
            end

            num_active_houses = sum(current_data.flow_kg_h >= absolute_flow_floor_kg_h);
            is_csac_active = (num_active_houses >= min_active_houses);

            [T_junction_ukf_C, ukf_master_offset] = calculate_main_pipe_temp( ...
                current_data, current_T_soil_C, U_csac, absolute_flow_floor_kg_h, ukf_states);
            [T_junction_pf_C, pf_master_offset] = calculate_main_pipe_temp( ...
                current_data, current_T_soil_C, U_csac, absolute_flow_floor_kg_h, pf_states);
            
            if all(isnan(T_junction_ukf_C)) || all(isnan(T_junction_pf_C))
                continue
            end

            current_data.T_main_ukf_C = T_junction_ukf_C;
            current_data.T_main_pf_C  = T_junction_pf_C;
            
            logger_ukf.timestamps(t) = time;
            logger_pf.timestamps(t)  = time;

            % ============================================================
            % PER-HOUSE LOOP
            % ============================================================
            for i = 1:num_houses_csac
                log_ukf = false;
                log_pf  = false;

                house_id = houses_csac_ids(i);
                house_data = current_data(current_data.house_id == house_id, :);
                if isempty(house_data)
                    continue
                end

                % --- Flow validity gates ---
                [is_flow_sufficient_ukf, alpha_ukf, ~, ~] = flow_validity_gate( ...
                    house_data.flow_kg_h, ukf_states{i}.x(2), house_data.length_service_m, alpha_min);
                is_flow_sufficient_ukf = is_flow_sufficient_ukf && (house_data.flow_kg_h > absolute_flow_floor_kg_h);

                [is_flow_sufficient_pf, alpha_pf, ~, ~] = flow_validity_gate( ...
                    house_data.flow_kg_h, pf_states{i}.x(2), house_data.length_service_m, alpha_min);
                is_flow_sufficient_pf = is_flow_sufficient_pf && (house_data.flow_kg_h > absolute_flow_floor_kg_h);

                % --- Delta-T gates ---
                delta_T_ukf_sufficient = ((house_data.T_main_ukf_C - current_T_soil_C) >= delta_T_gate_threshold);
                delta_T_pf_sufficient  = ((house_data.T_main_pf_C  - current_T_soil_C) >= delta_T_gate_threshold);

                % --- Combined pre-condition ---
                can_update_ukf = is_csac_active && is_flow_sufficient_ukf && delta_T_ukf_sufficient;
                can_update_pf  = is_csac_active && is_flow_sufficient_pf  && delta_T_pf_sufficient;

                if csac == 1 && t <= 20
                    print_full_gate_summary(csac, time, house_id, ...
                        is_csac_active, is_flow_sufficient_ukf, delta_T_ukf_sufficient, ...
                        house_data.T_main_ukf_C, current_T_soil_C, house_data.flow_kg_h, ...
                        alpha_ukf, 'UKF');
                end

                % ========================================================
                % UKF UPDATE
                % ========================================================
                if can_update_ukf

                    % --- Hibernation reset ---
                    hours_since_last_ukf = Inf;
                    if ~isnat(last_update_timestamp_ukf(i))
                        hours_since_last_ukf = hours(time - last_update_timestamp_ukf(i));
                    end

                    if hours_since_last_ukf > hibernation_reset_threshold_hours
                        last_valid_T_main_ukf_C(i) = NaN;
                        ukf_states{i}.P = ukf_states{i}.P * 2.0;
                        ukf_innovation_gate(i) = innovation_gate_initial;
                        if debug_print_update_result
                            fprintf('UKF HIBERNATION RESET | house %d | time %s | %.0f hrs since last update\n', ...
                                house_id, string(time), hours_since_last_ukf);
                        end
                    end

                    % --- Stability gate ---
                    if isnan(last_valid_T_main_ukf_C(i))
                        delta_T_main_change_ukf = 0;
                    else
                        delta_T_main_change_ukf = abs(house_data.T_main_ukf_C - last_valid_T_main_ukf_C(i));
                    end

                    if debug_disable_stability_gate
                        is_system_stable = true;
                    else
                        is_system_stable = (delta_T_main_change_ukf < max_delta_T_change);
                    end

                    if is_system_stable
                        % --- Compute innovation ---
                        predicted_temp = get_supply_temp(house_data.T_main_ukf_C, ...
                            house_data.flow_kg_h, ukf_states{i}.x(2), ...
                            house_data.length_service_m, current_T_soil_C);
                        innovation = house_data.T_supply_C - (predicted_temp - ukf_states{i}.x(1));

                        % --- DUAL INNOVATION GATE ---
                        effective_ukf_gate = max(ukf_innovation_gate(i), 3.0);
                        passes_adaptive_gate   = abs(innovation) < effective_ukf_gate;
                        passes_hard_gate       = abs(innovation) < max_innovation;
                        passes_innovation_gate = passes_adaptive_gate && passes_hard_gate;

                        if passes_innovation_gate
                            x_before = ukf_states{i}.x;
                            P_before = ukf_states{i}.P;

                            [ukf_states{i}, diagnostics_ukf] = update_ukf_house( ...
                                ukf_states{i}, house_data, current_T_soil_C, config);

                            if any(isnan(ukf_states{i}.x), 'all') || any(isnan(ukf_states{i}.P), 'all')
                                warning('UKF returned NaN for house %d at %s. Restoring previous state.', ...
                                    house_id, string(time));
                                ukf_states{i}.x = x_before;
                                ukf_states{i}.P = P_before;
                            else
                                ukf_offsets(i) = ukf_states{i}.x(1);
                                log_ukf = true;
                                last_valid_T_main_ukf_C(i) = house_data.T_main_ukf_C;
                                last_update_timestamp_ukf(i) = time;
                            end

                            % Update adaptive gate width
                            if isfield(diagnostics_ukf, 'P_zz') && isfinite(diagnostics_ukf.P_zz) && diagnostics_ukf.P_zz > 0
                                raw_gate = sqrt(diagnostics_ukf.P_zz) * gate_N_sigma;
                                ukf_innovation_gate(i) = max(1.0, min(raw_gate, 5.0));  % Cap at 5°C
                            else
                                warning('Invalid UKF P_zz for house %d at %s', house_id, string(time));
                            end

                            if debug_print_update_result
                                fprintf('UKF UPDATE | house %d | time %s | innov=%.3f | dx=[%.4f %.4f]\n', ...
                                    house_id, string(time), innovation, ...
                                    ukf_states{i}.x(1) - x_before(1), ukf_states{i}.x(2) - x_before(2));
                            end

                        else
                            % Gate rejection logging
                            if debug_print_update_result
                                fprintf('UKF GATE REJECT | house %d | time %s | innov=%.3f | adaptive=%.3f | hard=%.3f\n', ...
                                    house_id, string(time), innovation, effective_ukf_gate, max_innovation);
                            end
                        end % innovation gate

                    else
                        % Stability gate rejection logging
                        if debug_print_update_result
                            fprintf('UKF STABILITY REJECT | house %d | time %s | dT_main=%.3f > threshold=%.3f\n', ...
                                house_id, string(time), delta_T_main_change_ukf, max_delta_T_change);
                        end
                    end % stability gate

                end % can_update_ukf

                % ========================================================
                % PF UPDATE
                % ========================================================
                if can_update_pf

                    % --- Hibernation reset ---
                    hours_since_last_pf = Inf;
                    if ~isnat(last_update_timestamp_pf(i))
                        hours_since_last_pf = hours(time - last_update_timestamp_pf(i));
                    end

                    if hours_since_last_pf > hibernation_reset_threshold_hours
                        last_valid_T_main_pf_C(i) = NaN;
                        pf_particles{i} = initialize_pf_state(num_particles, pf_states{i}.x, P_base, config);
                        pf_innovation_gate(i) = innovation_gate_initial;
                        if debug_print_update_result
                            fprintf('PF HIBERNATION RESET | house %d | time %s | %.0f hrs since last update\n', ...
                                house_id, string(time), hours_since_last_pf);
                        end
                    end

                    % --- Stability gate ---
                    if isnan(last_valid_T_main_pf_C(i))
                        delta_T_main_change_pf = 0;
                    else
                        delta_T_main_change_pf = abs(house_data.T_main_pf_C - last_valid_T_main_pf_C(i));
                    end

                    if debug_disable_stability_gate
                        is_system_stable = true;
                    else
                        is_system_stable = (delta_T_main_change_pf < max_delta_T_change);
                    end

                    if is_system_stable
                        predicted_temp = get_supply_temp(house_data.T_main_pf_C, ...
                            house_data.flow_kg_h, pf_states{i}.x(2), ...
                            house_data.length_service_m, current_T_soil_C);
                        innovation = house_data.T_supply_C - (predicted_temp - pf_states{i}.x(1));

                        % --- DUAL INNOVATION GATE ---
                        effective_pf_gate = max(pf_innovation_gate(i), 3.0);
                        passes_adaptive_gate   = abs(innovation) < effective_pf_gate;
                        passes_hard_gate       = abs(innovation) < max_innovation;
                        passes_innovation_gate = passes_adaptive_gate && passes_hard_gate;

                        if passes_innovation_gate
                            x_before = pf_states{i}.x;

                            [pf_particles{i}, est_pf, cov_pf, diagnostics_pf] = update_pf_house( ...
                                pf_particles{i}, house_data, current_T_soil_C, R_base, Q_base, config);

                            if any(isnan(est_pf), 'all') || any(isnan(cov_pf), 'all')
                                warning('PF returned NaN for house %d at %s', house_id, string(time));
                            else
                                pf_states{i}.x = est_pf;
                                pf_states{i}.P = cov_pf;
                                pf_offsets(i) = pf_states{i}.x(1);
                                log_pf = true;
                                last_valid_T_main_pf_C(i) = house_data.T_main_pf_C;
                                last_update_timestamp_pf(i) = time;

                                % Update adaptive gate width (PF)
                                if isfield(diagnostics_pf, 'P_zz') && isfinite(diagnostics_pf.P_zz) && diagnostics_pf.P_zz > 0
                                    raw_gate = sqrt(diagnostics_pf.P_zz) * gate_N_sigma;
                                    pf_innovation_gate(i) = max(1.0, min(raw_gate, 5.0));  % Cap at 5°C
                                else
                                    warning('Invalid PF P_zz for house %d at %s', house_id, string(time));
                                end

                                if debug_print_update_result
                                    fprintf('PF UPDATE | house %d | time %s | innov=%.3f | dx=[%.4f %.4f]\n', ...
                                        house_id, string(time), innovation, ...
                                        pf_states{i}.x(1) - x_before(1), pf_states{i}.x(2) - x_before(2));
                                end
                            end

                        else
                            if debug_print_update_result
                                fprintf('PF GATE REJECT | house %d | time %s | innov=%.3f | adaptive=%.3f | hard=%.3f\n', ...
                                    house_id, string(time), innovation, effective_pf_gate, max_innovation);
                            end
                        end % innovation gate

                    else
                        if debug_print_update_result
                            fprintf('PF STABILITY REJECT | house %d | time %s | dT_main=%.3f > threshold=%.3f\n', ...
                                house_id, string(time), delta_T_main_change_pf, max_delta_T_change);
                        end
                    end % stability gate

                end % can_update_pf

                % --- Logging ---
                if log_ukf
                    logger_ukf = update_logger(logger_ukf, t, i, time, ukf_states{i}, diagnostics_ukf);
                elseif t > 1
                    logger_ukf.state_estimates(:, i, t) = logger_ukf.state_estimates(:, i, t-1);
                    logger_ukf.covariance_posterior(:, i, t) = logger_ukf.covariance_posterior(:, i, t-1);
                end
                if log_pf
                    logger_pf = update_logger(logger_pf, t, i, time, pf_states{i}, diagnostics_pf);
                elseif t > 1
                    logger_pf.state_estimates(:, i, t) = logger_pf.state_estimates(:, i, t-1);
                    logger_pf.covariance_posterior(:, i, t) = logger_pf.covariance_posterior(:, i, t-1);
                end

            end % per-house loop

            % ============================================================
            % POST-HOUSE-LOOP: Apply damped master_offset corrections
            % ============================================================
            if ~isnan(ukf_master_offset) && abs(ukf_master_offset) > master_offset_threshold
                correction_ukf = master_offset_lr * ukf_master_offset;
                for j = 1:num_houses_csac
                    ukf_states{j}.x(1) = ukf_states{j}.x(1) - correction_ukf;
                end
                if debug_print_update_result
                    fprintf('UKF MASTER OFFSET | time %s | raw=%.4f | applied=%.4f\n', ...
                        string(time), ukf_master_offset, correction_ukf);
                end
            end
            
            if ~isnan(pf_master_offset) && abs(pf_master_offset) > master_offset_threshold
                correction_pf = master_offset_lr * pf_master_offset;
                for j = 1:num_houses_csac
                    pf_states{j}.x(1) = pf_states{j}.x(1) - correction_pf;
                    pf_particles{j}(1,:) = pf_particles{j}(1,:) - correction_pf;
                end
                if debug_print_update_result
                    fprintf('PF MASTER OFFSET | time %s | raw=%.4f | applied=%.4f\n', ...
                        string(time), pf_master_offset, correction_pf);
                end
            end

            % --- Mean centering (currently disabled) ---
            if error_meaning && ~debug_disable_mean_centering
                mean_ukf_offsets = mean(ukf_offsets, 'omitmissing');
                mean_pf_offsets  = mean(pf_offsets,  'omitmissing');
                for j = 1:num_houses_csac
                    if (~isnan(mean_ukf_offsets) && ~isnan(ukf_states{j}.x(1)))
                        ukf_states{j}.x(1) = ukf_states{j}.x(1) - mean_ukf_offsets;
                    end
                    if (~isnan(mean_pf_offsets) && ~isnan(pf_states{j}.x(1)))
                        pf_states{j}.x(1) = pf_states{j}.x(1) - mean_pf_offsets;
                        pf_particles{j}(1,:) = pf_particles{j}(1,:) - mean_pf_offsets;
                    end
                end
            end

        end % time loop

        plot_diagnostics(logger_ukf, ground_truth_csac, csac, network, output_folder_ukf);
        plot_diagnostics(logger_pf,  ground_truth_csac, csac, network, output_folder_pf);
        save_logger_to_csv(logger_ukf, output_folder_ukf, strcat('ukf_csac_', string(csac)));
        save_logger_to_csv(logger_pf,  output_folder_pf,  strcat('pf_csac_', string(csac)));
    end % csac loop
end % network loop

close(w);
disp('Analysis complete. All diagnostic plots have been saved.')