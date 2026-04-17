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

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'src/gates', 'config')

% read config-file
config = jsondecode(fileread("config.json"));

% --- DEBUG SETTINGS (from config) ---
debug_disable_stability_gate = config.project.debug.disable_stability_gate;
debug_disable_mean_centering = config.project.debug.disable_mean_centering;
debug_print_update_result = config.project.debug.print_update_result;

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
        ukf_innovation_gate = nan([1,num_houses_csac]);
        
        last_valid_T_main_ukf_C = nan(num_houses_csac,1);
        last_update_timestamp_ukf = NaT(num_houses_csac,1);

        % Gate statistics
        gate_reject_count_ukf = 0;
        gate_accept_count_ukf = 0;

        for i = 1:num_houses_csac
            x_init = [randn()*0.3; 0.12 + randn()*0.03];
            ukf_states{i} = initialize_kalman_filter(x_init, Q_base, R_base, P_base);
            ukf_innovation_gate(i) = innovation_gate_initial;
        end
        
        % set up loggers
        logger_ukf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);

        % === LOG INITIAL STATES ===
        initial_states_ukf = zeros(2, num_houses_csac);
        for i = 1:num_houses_csac
            initial_states_ukf(:, i) = ukf_states{i}.x;
            logger_ukf.state_estimates(:, i, 1) = ukf_states{i}.x;
            logger_ukf.covariance_posterior(:, i, 1) = [ukf_states{i}.P(1,1); ukf_states{i}.P(2,2)];
        end

        ukf_offsets = nan([num_houses_csac,1]);

        % === HEATING SEASON STATE ===
        season_state = initialize_season_state();
        
        % Snapshot of last known good states
        ukf_states_snapshot = cell(size(ukf_states));
        for i = 1:num_houses_csac
            ukf_states_snapshot{i}.x = ukf_states{i}.x;
            ukf_states_snapshot{i}.P = ukf_states{i}.P;
        end

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
            %% HEATING SEASON GATE
            %% ============================================================
            current_date = dateshift(time, 'start', 'day');
            [season_state, actions] = manage_heating_season(...
                season_state, current_date, t, daily_T_air_max_table, config, csac);

            % Apply snapshot restore/save
            if actions.do_restore_snapshot || actions.do_save_snapshot
                [ukf_states, ukf_states_snapshot] = apply_season_actions(...
                    actions, ukf_states, ukf_states_snapshot, P_base, config);
            end

            % Track last active date for gap calculation
            if actions.season_is_active
                season_state.last_active_date = current_date;
            end

            % Skip if not in season
            if actions.skip_timestep
                if t > 1
                    for i = 1:num_houses_csac
                        logger_ukf.state_estimates(:, i, t) = logger_ukf.state_estimates(:, i, t-1);
                        logger_ukf.covariance_posterior(:, i, t) = logger_ukf.covariance_posterior(:, i, t-1);
                    end
                end
                logger_ukf.timestamps(t) = time;
                continue;
            end

            %% ============================================================
            %% ACTIVE SEASON: Normal processing continues below
            %% ============================================================

            num_active_houses = sum(current_data.flow_kg_h >= absolute_flow_floor_kg_h);
            is_csac_active = (num_active_houses >= config.project.initialization.min_active_houses);

            [T_junction_ukf_C, ukf_master_offset] = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, absolute_flow_floor_kg_h, ukf_states);
            
            if all(isnan(T_junction_ukf_C))
                continue
            end

            current_data.T_main_ukf_C = T_junction_ukf_C;
            
            logger_ukf.timestamps(t) = time;

            for i = 1:num_houses_csac
                log_ukf = false;
                house_id = houses_csac_ids(i);
                house_data = current_data(current_data.house_id==house_id,:);
                if isempty(house_data)
                    continue
                end
                [is_flow_sufficient_ukf, alpha_ukf, theta_ukf, min_flow_ukf_kg_h] = flow_validity_gate(house_data.flow_kg_h, ukf_states{i}.x(2), house_data.length_service_m, alpha_min);
                is_flow_sufficient_ukf = is_flow_sufficient_ukf && (house_data.flow_kg_h>absolute_flow_floor_kg_h);

                delta_T_ukf_sufficient = ((house_data.T_main_ukf_C - current_T_soil_C) >= delta_T_gate_threshold);

                can_update_ukf = is_csac_active && is_flow_sufficient_ukf && delta_T_ukf_sufficient;

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

                if log_ukf
                    logger_ukf = update_logger(logger_ukf, t, i, time, ukf_states{i}, diagnostics_ukf);
                elseif t > 1
                    logger_ukf.state_estimates(:, i, t) = logger_ukf.state_estimates(:, i, t-1);
                    logger_ukf.covariance_posterior(:, i, t) = logger_ukf.covariance_posterior(:, i, t-1);
                end
            end

            %% ============================================================
            %% SNAPSHOT UPDATE (only after season has run long enough)
            %% ============================================================
            min_season_for_cooldown = config.project.heating_season_gate.min_season_days_for_cooldown;
            if ~isnat(season_state.start_date)
                current_season_days = days(current_date - season_state.start_date);
            else
                current_season_days = 0;
            end
            
            if current_season_days >= min_season_for_cooldown
                % Season is established — safe to update snapshot
                season_state.snapshot_timestep = t;
                for i = 1:num_houses_csac
                    ukf_states_snapshot{i}.x = ukf_states{i}.x;
                    ukf_states_snapshot{i}.P = ukf_states{i}.P;
                end
            end

            %% ============================================================
            %% MASTER OFFSET APPLICATION
            %% ============================================================
            if ~isnan(ukf_master_offset) && abs(ukf_master_offset) > master_offset_deadzone
                if num_active_houses >= min_active_houses
                    max_master_offset = 0.5;
                    clamped_offset = max(-max_master_offset, min(max_master_offset, ukf_master_offset));
                    damped_ukf_offset = master_offset_gamma * clamped_offset;
                    for i = 1:num_houses_csac
                        ukf_states{i}.x(1) = ukf_states{i}.x(1) - damped_ukf_offset;
                    end
                end
            end
        end

        % Print gate and season statistics
        fprintf('\n=== CSAC %d Final Statistics ===\n', csac);
        fprintf('UKF: %d accepted, %d rejected (%.1f%% rejection rate)\n', ...
            gate_accept_count_ukf, gate_reject_count_ukf, ...
            100*gate_reject_count_ukf / max(1, gate_accept_count_ukf + gate_reject_count_ukf));
        fprintf('Season gate: %d active days, %d inactive days (%.1f%% active)\n', ...
            season_state.active_days, season_state.inactive_days, ...
            100*season_state.active_days / max(1, season_state.active_days + season_state.inactive_days));
        fprintf('Heating seasons started: %d\n', season_state.count);
        if ~isnan(season_state.end_timestep)
            fprintf('Last season ended at timestep %d\n', season_state.end_timestep);
        end
        if season_state.cooldown_active
            fprintf('Cooldown active from real season ending %s\n', string(season_state.last_real_season_end_date));
        end
        fprintf('=====================================\n\n');

        plot_diagnostics(logger_ukf, ground_truth_csac, csac, network, output_folder_ukf);
        save_logger_to_csv(logger_ukf, output_folder_ukf, strcat('ukf_csac_', string(csac)));

        save_diagnostic_summary(logger_ukf, ground_truth_csac, csac, output_folder_ukf, 'ukf');
        save_diagnostic_summary_detailed(logger_ukf, ground_truth_csac, csac, output_folder_ukf, 'ukf');
        save_daily_diagnostics(logger_ukf, ground_truth_csac, csac, output_folder_ukf, 'ukf');
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
    % Waitbar already closed or invalid
end
disp('Analysis complete. All diagnostic plots have been saved.')