% main.m
% Clear workspace, add paths, set parameters
clear all
close all
clc

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'src/gates', 'config')

% read config-file
config = jsondecode(fileread("config.json"));

% --- DEBUG SETTINGS (from config) ---
debug_disable_stability_gate = config.project.debug.disable_stability_gate;
debug_disable_mean_centering = config.project.debug.disable_mean_centering;
debug_print_update_result = config.project.debug.print_update_result;

% --- HEATING SEASON SETTINGS (from config) ---
season_nis_threshold = config.project.heating_season.nis_threshold;
season_window_days = config.project.heating_season.window_days;
season_consec_bad_to_end = config.project.heating_season.consec_bad_to_end;
season_min_active_fraction = config.project.heating_season.min_active_fraction;
season_restart_nis_threshold = config.project.heating_season.restart_nis_threshold;
season_consec_good_to_restart = config.project.heating_season.consec_good_to_restart;
season_P_inflation_factor = config.project.heating_season.P_inflation_factor;
season_min_update_rate = config.project.heating_season.min_update_rate;
season_warmup_days = config.project.heating_season.warmup_days;

% --- HEATING SEASON GATE SETTINGS (from config) ---
season_gate_T_air_max_threshold = config.project.heating_season_gate.T_air_max_threshold;
season_gate_lookback_days = config.project.heating_season_gate.lookback_days;
season_gate_P_restart_fraction = config.project.heating_season_gate.P_restart_fraction;

% --- MASTER OFFSET SETTINGS (from config) ---
master_offset_gamma = config.project.master_offset.gamma;
master_offset_deadzone = config.project.master_offset.deadzone;

% --- HARD INNOVATION GATE (from config) ---
hard_innovation_gate_C = config.project.initialization.max_innovation_C;
max_nis = config.project.initialization.max_nis;
fprintf('Hard innovation gate set to %.2f °C\n', hard_innovation_gate_C);

% ukf-configuration
R_base = config.project.initialization.ukf.measurement_noise^2;
Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, (config.project.initialization.ukf.process_noise_U)^2]);
P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, (config.project.initialization.ukf.state_uncertainty_U)^2]);
innovation_gate_initial = config.project.initialization.ukf.state_uncertainty_offset * config.project.initialization.innovation_gate_N_sigma;

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

networks = config.project.datasets.datasets;

% load network topology and measurement data
[T_soil_C, T_air_C] = soilTemp(config);

%% ============================================================
%% PRE-COMPUTE DAILY T_AIR_MAX LOOKUP TABLE
%% ============================================================
fprintf('Pre-computing daily T_air_max lookup table...\n');

% Extract unique dates from air temperature data
T_air_C.date = dateshift(T_air_C.time, 'start', 'day');
[air_date_groups, air_unique_dates] = findgroups(T_air_C.date);

% Compute daily max air temperature
daily_T_air_max_values = splitapply(@max, T_air_C.values, air_date_groups);
daily_T_air_max_table = table(air_unique_dates, daily_T_air_max_values, ...
    'VariableNames', {'date', 'T_air_max'});
daily_T_air_max_table = sortrows(daily_T_air_max_table, 'date');

fprintf('  Built daily T_air_max table: %d days (%.1f to %.1f °C)\n', ...
    height(daily_T_air_max_table), ...
    min(daily_T_air_max_table.T_air_max), ...
    max(daily_T_air_max_table.T_air_max));



output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'),'/ukf');
output_folder_pf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'),'/pf');
w = waitbar(0.0, "Starting analysis");

% run through all the networks from the config-file
for network = networks
    % load relevant data for the network
    [meter_data, network_data, topology] = importData(config,network);
    timestamps = unique(meter_data.timestamp);
    
    % find the ids of the csacs from the topology
    csac_ids = [topology.cul_de_sacs.id];
    prev_houses = 0;

    % loop through all cul-de-sacs in the network
    for csac=csac_ids
        
        % extract relevant data for the csac
        U_csac = topology.pipe_parameters.csac_pipe.insulation_W_m_K;
        houses_csac = ([topology.houses.cul_de_sac_id]==csac);
        topology_csac = topology.houses(houses_csac);
        num_houses_csac = sum(houses_csac);
        houses_csac_ids = sort([topology.houses(houses_csac).id]);
        true_offset = zeros([size(topology_csac),1]);

        ground_truth_csac = table([topology_csac.id]', true_offset, [topology_csac.service_pipe_insulation_W_m_K]', ...
            'VariableNames', {'house_id', 'true_offset', 'true_U'});
        
        csac_table = table([topology_csac.id]', [topology_csac.service_pipe_len_m]', [topology_csac.dist_on_cul_de_sac_m]', 'VariableNames',{'house_id','length_service_m','x_pos_m'});

        meter_data_csac = meter_data(ismember(meter_data.house_id,houses_csac_ids),:);
        meter_data_csac = join(meter_data_csac, csac_table, "Keys", "house_id");

        % initialize filters for each house
        ukf_states = cell(1, num_houses_csac);
        pf_particles = cell(1, num_houses_csac);
        pf_states = cell(1, num_houses_csac);
        ukf_innovation_gate = nan([1,num_houses_csac]);
        pf_innovation_gate = nan([1,num_houses_csac]);
        
        last_valid_T_main_ukf_C = nan(num_houses_csac,1);
        last_valid_T_main_pf_C = nan(num_houses_csac,1);
        last_update_timestamp_ukf = NaT(num_houses_csac,1);
        last_update_timestamp_pf = NaT(num_houses_csac,1);

        % Gate statistics
        gate_reject_count_ukf = 0;
        gate_reject_count_pf = 0;
        gate_accept_count_ukf = 0;
        gate_accept_count_pf = 0;

        for i = 1:num_houses_csac
            x_init = [randn()*0.3; 0.12 + randn()*0.03];
            ukf_states{i} = initialize_kalman_filter(x_init, Q_base, R_base, P_base);
            pf_particles{i} = initialize_pf_state(num_particles, x_init, P_base, config);
            pf_states{i}.x = x_init;
            pf_states{i}.P = P_base;
            ukf_innovation_gate(i) = innovation_gate_initial;
            pf_innovation_gate(i) = innovation_gate_initial;
        end
        
        % set up loggers
        logger_ukf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);
        logger_pf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);

        ukf_offsets = nan([num_houses_csac,1]);
        pf_offsets = nan([num_houses_csac,1]);

        % === HEATING SEASON STATE TRACKING ===
        season_active = false;  % Start inactive — wait for gate to confirm
        was_season_active = false;
        season_end_timestep = NaN;
        season_start_timestep = NaN;
        season_count = 0;
        
        % Snapshot of last known good states
        ukf_states_snapshot = cell(size(ukf_states));
        pf_states_snapshot = cell(size(pf_states));
        snapshot_timestep = NaN;
        for i = 1:num_houses_csac
            ukf_states_snapshot{i}.x = ukf_states{i}.x;
            ukf_states_snapshot{i}.P = ukf_states{i}.P;
            pf_states_snapshot{i}.x = pf_states{i}.x;
            pf_states_snapshot{i}.P = pf_states{i}.P;
        end

        % Season gate statistics
        season_gate_active_days = 0;
        season_gate_inactive_days = 0;
        last_gate_check_date = NaT;

        % run through the timesteps
        waitbar((csac)/length(csac_ids),w,strcat("Running network. At timestep: 0/", string(length(timestamps)),", csac: ",string(csac),"/",string(length(csac_ids)-1),", network: ", string(network)))
        for t = 1:length(timestamps)

            if mod(t,10)==0
                waitbar(t/(length(timestamps)*length(csac_ids))+(csac)/length(csac_ids),w,strcat("Running network. At timestep: ", string(t), "/", string(length(timestamps)),", csac: ",string(csac),"/",string(length(csac_ids)-1),", network: ", string(network)))
            end
            time = timestamps(t);

            current_data = meter_data_csac(meter_data_csac.timestamp==time,:);
            current_data = sortrows(current_data, 'house_id');

            current_T_soil_C = T_soil_C(T_soil_C.time==time,:).values;
            current_T_air_C = T_air_C((year(T_air_C.time)==year(time)) & (month(T_air_C.time)==month(time)) & (day(T_air_C.time)==day(time)),:).values;

            if isempty(current_data) || isempty(current_T_soil_C)
                continue
            end

            %% ============================================================
            %% HEATING SEASON GATE (temperature-based, checked once per day)
            %% ============================================================
            current_date = dateshift(time, 'start', 'day');
            
            if isnat(last_gate_check_date) || current_date > last_gate_check_date
                % New day — evaluate the gate
                last_gate_check_date = current_date;
                
                [gate_is_active, recent_temps, n_valid_days] = check_heating_season_gate(...
                    current_date, season_gate_lookback_days, ...
                    season_gate_T_air_max_threshold, daily_T_air_max_table);
                
                % --- Detect transitions ---
                was_season_active = season_active;
                
                if gate_is_active && ~was_season_active
                    % === SEASON START ===
                    season_active = true;
                    season_count = season_count + 1;
                    season_start_timestep = t;
                    
                    fprintf('\n*** CSAC %d: HEATING SEASON #%d START at %s (timestep %d) ***\n', ...
                        csac, season_count, string(time), t);
                    fprintf('    Last %d days T_air_max: [%s] (all < %.0f°C)\n', ...
                        season_gate_lookback_days, ...
                        strjoin(arrayfun(@(x) sprintf('%.1f', x), recent_temps(isfinite(recent_temps)), 'UniformOutput', false), ', '), ...
                        season_gate_T_air_max_threshold);
                    
                    if season_count == 1
                        % First season: use initial P (cold start)
                        fprintf('    First season — using initial P (cold start)\n');
                    else
                        % Subsequent season: inflate P from snapshot
                        fprintf('    Restarting from snapshot (timestep %d)\n', snapshot_timestep);
                        for i = 1:num_houses_csac
                            % Restore point estimates from snapshot
                            ukf_states{i}.x = ukf_states_snapshot{i}.x;
                            pf_states{i}.x = pf_states_snapshot{i}.x;
                            
                            % Inflate covariance: blend between snapshot P and P_base
                            P_snap = ukf_states_snapshot{i}.P;
                            P_restart = P_snap + season_gate_P_restart_fraction * (P_base - P_snap);
                            % Ensure we don't exceed P_base
                            P_restart(1,1) = min(P_restart(1,1), P_base(1,1));
                            P_restart(2,2) = min(P_restart(2,2), P_base(2,2));
                            % Ensure we don't go below snapshot
                            P_restart(1,1) = max(P_restart(1,1), P_snap(1,1));
                            P_restart(2,2) = max(P_restart(2,2), P_snap(2,2));
                            
                            ukf_states{i}.P = P_restart;
                            pf_states{i}.P = P_restart;
                            
                            if i == 1
                                fprintf('    House %d: P_snap=[%.4f, %.6f], P_restart=[%.4f, %.6f], P_base=[%.4f, %.6f]\n', ...
                                    houses_csac_ids(i), ...
                                    sqrt(P_snap(1,1)), sqrt(P_snap(2,2)), ...
                                    sqrt(P_restart(1,1)), sqrt(P_restart(2,2)), ...
                                    sqrt(P_base(1,1)), sqrt(P_base(2,2)));
                            end
                        end
                    end
                    
                elseif ~gate_is_active && was_season_active
                    % === SEASON END (gate-based) ===
                    season_active = false;
                    season_end_timestep = t;
                    
                    fprintf('\n*** CSAC %d: HEATING SEASON #%d END at %s (timestep %d) ***\n', ...
                        csac, season_count, string(time), t);
                    fprintf('    T_air_max exceeded %.0f°C in lookback window: [%s]\n', ...
                        season_gate_T_air_max_threshold, ...
                        strjoin(arrayfun(@(x) sprintf('%.1f', x), recent_temps(isfinite(recent_temps)), 'UniformOutput', false), ', '));
                    
                    % Save snapshot of current states
                    for i = 1:num_houses_csac
                        ukf_states_snapshot{i}.x = ukf_states{i}.x;
                        ukf_states_snapshot{i}.P = ukf_states{i}.P;
                        pf_states_snapshot{i}.x = pf_states{i}.x;
                        pf_states_snapshot{i}.P = pf_states{i}.P;
                    end
                    snapshot_timestep = t;
                    fprintf('    Snapshot saved at timestep %d\n', t);
                end
                
                % Track gate statistics
                if season_active
                    season_gate_active_days = season_gate_active_days + 1;
                else
                    season_gate_inactive_days = season_gate_inactive_days + 1;
                end
            end
            
            % --- If season is not active, carry forward and skip ---
            if ~season_active
                if t > 1
                    for i = 1:num_houses_csac
                        logger_ukf.state_estimates(:, i, t) = logger_ukf.state_estimates(:, i, t-1);
                        logger_ukf.covariance_posterior(:, i, t) = logger_ukf.covariance_posterior(:, i, t-1);
                        logger_pf.state_estimates(:, i, t) = logger_pf.state_estimates(:, i, t-1);
                        logger_pf.covariance_posterior(:, i, t) = logger_pf.covariance_posterior(:, i, t-1);
                    end
                end
                logger_ukf.timestamps(t) = time;
                logger_pf.timestamps(t) = time;
                continue;  % Skip to next timestep
            end

            %% ============================================================
            %% ACTIVE SEASON: Normal processing continues below
            %% ============================================================

            num_active_houses = sum(current_data.flow_kg_h >= absolute_flow_floor_kg_h);
            is_csac_active = (num_active_houses >= config.project.initialization.min_active_houses);

            [T_junction_ukf_C, ukf_master_offset] = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, absolute_flow_floor_kg_h, ukf_states);
            [T_junction_pf_C, pf_master_offset] = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, absolute_flow_floor_kg_h, pf_states);
            
            if all(isnan(T_junction_ukf_C)) || all(isnan(T_junction_pf_C))
                continue
            end

            current_data.T_main_ukf_C = T_junction_ukf_C;
            current_data.T_main_pf_C = T_junction_pf_C;
            
            logger_ukf.timestamps(t) = time;
            logger_pf.timestamps(t) = time;

            for i = 1:num_houses_csac
                log_ukf = false;
                log_pf = false;
                house_id = houses_csac_ids(i);
                house_data = current_data(current_data.house_id==house_id,:);
                if isempty(house_data)
                    continue
                end
                [is_flow_sufficient_ukf, alpha_ukf, theta_ukf, min_flow_ukf_kg_h] = flow_validity_gate(house_data.flow_kg_h, ukf_states{i}.x(2), house_data.length_service_m, alpha_min);
                is_flow_sufficient_ukf = is_flow_sufficient_ukf && (house_data.flow_kg_h>absolute_flow_floor_kg_h);

                [is_flow_sufficient_pf, alpha_pf, theta_pf, min_flow_pf_kg_h] = flow_validity_gate(house_data.flow_kg_h, pf_states{i}.x(2), house_data.length_service_m, alpha_min);
                is_flow_sufficient_pf = is_flow_sufficient_pf && (house_data.flow_kg_h>absolute_flow_floor_kg_h);
               
                delta_T_ukf_sufficient = ((house_data.T_main_ukf_C - current_T_soil_C) >= delta_T_gate_threshold);
                delta_T_pf_sufficient = ((house_data.T_main_pf_C - current_T_soil_C) >= delta_T_gate_threshold);

                can_update_ukf = is_csac_active && is_flow_sufficient_ukf && delta_T_ukf_sufficient;
                can_update_pf = is_csac_active && is_flow_sufficient_pf && delta_T_pf_sufficient;

                %% ============================================================
                %% UKF UPDATE BLOCK
                %% ============================================================
                if can_update_ukf
                    if isnan(last_valid_T_main_ukf_C(i))
                        delta_T_main_change_ukf = 0;
                    else
                        delta_T_main_change_ukf = abs(house_data.T_main_ukf_C-last_valid_T_main_ukf_C(i));
                    end

                    if debug_disable_stability_gate
                        is_system_stable = true;
                    else
                        is_system_stable = (delta_T_main_change_ukf < max_delta_T_change);
                    end

                    if is_system_stable
                        predicted_temp = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, ukf_states{i}.x(2), house_data.length_service_m, current_T_soil_C);
                        innovation = house_data.T_supply_C - (predicted_temp - ukf_states{i}.x(1));

                        if abs(innovation) > hard_innovation_gate_C
                            gate_reject_count_ukf = gate_reject_count_ukf + 1;
                            if debug_print_update_result
                                fprintf('UKF GATE REJECT | house %d | time %s | innovation=%.3f°C > gate=%.1f°C\n', ...
                                    house_id, string(time), innovation, hard_innovation_gate_C);
                            end
                        else
                            x_before = ukf_states{i}.x;
                            P_before = ukf_states{i}.P;

                            [ukf_states{i}, diagnostics_ukf] = update_ukf_house(ukf_states{i}, house_data, current_T_soil_C, config);
                            if any(isnan(ukf_states{i}.x), 'all') || any(isnan(ukf_states{i}.P), 'all')
                                warning('UKF returned NaN for house %d at %s', house_id, string(time));
                            else
                                % --- NIS gate: revert if update quality is poor ---
                                nis_value = NaN;
                                if isfield(diagnostics_ukf, 'y') && isfield(diagnostics_ukf, 'P_zz') ...
                                        && isfinite(diagnostics_ukf.P_zz) && diagnostics_ukf.P_zz > 0
                                    nis_value = (diagnostics_ukf.y^2) / diagnostics_ukf.P_zz;
                                end
                                
                                if isfinite(nis_value) && nis_value > max_nis
                                    % Revert state
                                    ukf_states{i}.x = x_before;
                                    ukf_states{i}.P = P_before;
                                    gate_reject_count_ukf = gate_reject_count_ukf + 1;
                                    if debug_print_update_result
                                        fprintf('UKF NIS REJECT | house %d | time %s | NIS=%.2f > %.1f | innov=%.3f\n', ...
                                            house_id, string(time), nis_value, max_nis, diagnostics_ukf.y);
                                    end
                                else
                                    % Accept update
                                    gate_accept_count_ukf = gate_accept_count_ukf + 1;
                                    ukf_offsets(i) = ukf_states{i}.x(1);
                                    log_ukf = true;
                                    last_valid_T_main_ukf_C(i) = house_data.T_main_ukf_C;
                                    last_update_timestamp_ukf(i) = time;
                                    
                                    if isfield(diagnostics_ukf, 'P_zz') && isfinite(diagnostics_ukf.P_zz) && diagnostics_ukf.P_zz > 0
                                        ukf_innovation_gate(i) = max(1, sqrt(diagnostics_ukf.P_zz) * config.project.initialization.innovation_gate_N_sigma);
                                    end
                                    if debug_print_update_result
                                        fprintf('UKF UPDATE | house %d | time %s | NIS=%.2f | innov=%.3f | dx=[%.4f %.6f]\n', ...
                                            house_id, string(time), nis_value, diagnostics_ukf.y, ...
                                            ukf_states{i}.x(1)-x_before(1), ukf_states{i}.x(2)-x_before(2));
                                    end
                                end
                            end
                        end
                    end
                end
                
                %% ============================================================
                %% PF UPDATE BLOCK
                %% ============================================================
                if can_update_pf
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
                        predicted_temp = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, pf_states{i}.x(2), house_data.length_service_m, current_T_soil_C);
                        innovation = house_data.T_supply_C - (predicted_temp - pf_states{i}.x(1));

                        if abs(innovation) > hard_innovation_gate_C
                            gate_reject_count_pf = gate_reject_count_pf + 1;
                            if debug_print_update_result
                                fprintf('PF  GATE REJECT | house %d | time %s | innovation=%.3f°C > gate=%.1f°C\n', ...
                                    house_id, string(time), innovation, hard_innovation_gate_C);
                            end
                        else
                            x_before = pf_states{i}.x;

                            [pf_particles{i}, est_pf, cov_pf, diagnostics_pf] = update_pf_house(pf_particles{i}, house_data, current_T_soil_C, R_base, Q_base, config);
                            
                            if any(isnan(est_pf), 'all') || any(isnan(cov_pf), 'all')
                                warning('PF returned NaN for house %d at %s', house_id, string(time));
                            else
                                % --- NIS gate ---
                                nis_value = NaN;
                                if isfield(diagnostics_pf, 'P_zz') && isfinite(diagnostics_pf.P_zz) && diagnostics_pf.P_zz > 0
                                    if isfield(diagnostics_pf, 'y') && isfinite(diagnostics_pf.y)
                                        nis_value = (diagnostics_pf.y^2) / diagnostics_pf.P_zz;
                                    end
                                end
                                
                                if isfinite(nis_value) && nis_value > max_nis
                                    pf_states{i}.x = x_before;
                                    gate_reject_count_pf = gate_reject_count_pf + 1;
                                    if debug_print_update_result
                                        fprintf('PF  NIS REJECT | house %d | time %s | NIS=%.2f > %.1f | innov=%.3f\n', ...
                                            house_id, string(time), nis_value, max_nis, diagnostics_pf.y);
                                    end
                                else
                                    gate_accept_count_pf = gate_accept_count_pf + 1;
                                    pf_states{i}.x = est_pf;
                                    pf_states{i}.P = cov_pf;
                                    pf_offsets(i) = pf_states{i}.x(1);
                                    log_pf = true;
                                    last_valid_T_main_pf_C(i) = house_data.T_main_pf_C;

                                    if isfield(diagnostics_pf, 'P_zz') && isfinite(diagnostics_pf.P_zz) && diagnostics_pf.P_zz > 0
                                        pf_innovation_gate(i) = max(1, sqrt(diagnostics_pf.P_zz)*config.project.initialization.innovation_gate_N_sigma);
                                    end
                                    if debug_print_update_result
                                        fprintf('PF  UPDATE | house %d | time %s | NIS=%.2f | innov=%.3f | dx=[%.4f %.6f]\n', ...
                                            house_id, string(time), nis_value, diagnostics_pf.y, ...
                                            pf_states{i}.x(1)-x_before(1), pf_states{i}.x(2)-x_before(2));
                                    end
                                end
                            end
                        end
                    end
                end

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
            end

            %% ============================================================
            %% SNAPSHOT UPDATE (always update during active season)
            %% ============================================================
            snapshot_timestep = t;
            for i = 1:num_houses_csac
                ukf_states_snapshot{i}.x = ukf_states{i}.x;
                ukf_states_snapshot{i}.P = ukf_states{i}.P;
                pf_states_snapshot{i}.x = pf_states{i}.x;
                pf_states_snapshot{i}.P = pf_states{i}.P;
            end

            %% ============================================================
            %% MASTER OFFSET APPLICATION
            %% ============================================================
            if ~isnan(ukf_master_offset) && abs(ukf_master_offset) > master_offset_deadzone
                damped_ukf_offset = master_offset_gamma * ukf_master_offset;
                for i = 1:num_houses_csac
                    ukf_states{i}.x(1) = ukf_states{i}.x(1) - damped_ukf_offset;
                end
            end

            if ~isnan(pf_master_offset) && abs(pf_master_offset) > master_offset_deadzone
                damped_pf_offset = master_offset_gamma * pf_master_offset;
                for i = 1:num_houses_csac
                    pf_states{i}.x(1) = pf_states{i}.x(1) - damped_pf_offset;
                    pf_particles{i}(1,:) = pf_particles{i}(1,:) - damped_pf_offset;
                end
            end

            if error_meaning && ~debug_disable_mean_centering
                mean_ukf_offsets = mean(ukf_offsets, 'omitmissing');
                mean_pf_offsets = mean(pf_offsets, 'omitmissing');
                for i=1:num_houses_csac
                    if (~isnan(mean_ukf_offsets) && ~isnan(ukf_states{i}.x(1)))
                        ukf_states{i}.x(1)=ukf_states{i}.x(1)-mean_ukf_offsets;
                    end
                    if (~isnan(mean_pf_offsets) && ~isnan(pf_states{i}.x(1)))
                        pf_states{i}.x(1)=pf_states{i}.x(1)-mean_pf_offsets;
                        pf_particles{i}(1,:) = pf_particles{i}(1,:) - mean_pf_offsets;
                    end
                end
            end
        end

        % Print gate and season statistics
        fprintf('\n=== CSAC %d Final Statistics ===\n', csac);
        fprintf('UKF: %d accepted, %d rejected (%.1f%% rejection rate)\n', ...
            gate_accept_count_ukf, gate_reject_count_ukf, ...
            100*gate_reject_count_ukf / max(1, gate_accept_count_ukf + gate_reject_count_ukf));
        fprintf('PF:  %d accepted, %d rejected (%.1f%% rejection rate)\n', ...
            gate_accept_count_pf, gate_reject_count_pf, ...
            100*gate_reject_count_pf / max(1, gate_accept_count_pf + gate_reject_count_pf));
        fprintf('Season gate: %d active days, %d inactive days (%.1f%% active)\n', ...
            season_gate_active_days, season_gate_inactive_days, ...
            100*season_gate_active_days / max(1, season_gate_active_days + season_gate_inactive_days));
        fprintf('Heating seasons started: %d\n', season_count);
        if ~isnan(season_end_timestep)
            fprintf('Last season ended at timestep %d\n', season_end_timestep);
        end
        fprintf('=====================================\n\n');

        plot_diagnostics(logger_ukf, ground_truth_csac, csac, network, output_folder_ukf);
        plot_diagnostics(logger_pf, ground_truth_csac, csac, network, output_folder_pf);
        save_logger_to_csv(logger_ukf, output_folder_ukf, strcat('ukf_csac_', string(csac)));
        save_logger_to_csv(logger_pf, output_folder_pf, strcat('pf_csac_', string(csac)));

        save_diagnostic_summary(logger_ukf, ground_truth_csac, csac, output_folder_ukf, 'ukf');
        save_diagnostic_summary(logger_pf, ground_truth_csac, csac, output_folder_pf, 'pf');
        save_diagnostic_summary_detailed(logger_ukf, ground_truth_csac, csac, output_folder_ukf, 'ukf');
        save_diagnostic_summary_detailed(logger_pf, ground_truth_csac, csac, output_folder_pf, 'pf');
        save_daily_diagnostics(logger_ukf, ground_truth_csac, csac, output_folder_ukf, 'ukf');
        save_daily_diagnostics(logger_pf, ground_truth_csac, csac, output_folder_pf, 'pf');
    end
   
end

close(w);
disp('Analysis complete. All diagnostic plots have been saved.')

% --- Helper function: check if heating season is active ---
    function [is_active, recent_max_temps, n_valid] = check_heating_season_gate(...
            current_date, lookback_days, threshold, daily_T_air_max_tbl)
        
        recent_max_temps = nan(lookback_days, 1);
        for d = 0:(lookback_days - 1)
            check_date = current_date - days(d);
            row_idx = find(daily_T_air_max_tbl.date == check_date, 1);
            if ~isempty(row_idx)
                recent_max_temps(d + 1) = daily_T_air_max_tbl.T_air_max(row_idx);
            end
        end
        
        n_valid = sum(isfinite(recent_max_temps));
        
        if n_valid < lookback_days
            is_active = false;
            return;
        end
        
        % ALL days in the lookback window must be below threshold
        is_active = all(recent_max_temps(isfinite(recent_max_temps)) < threshold);
    end