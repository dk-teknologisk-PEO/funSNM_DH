function cs = process_csac_timestep(cs, t, time, current_data, current_T_soil_C, ...
    daily_T_air_max_table, P_base, config, csac_id)
%PROCESS_CSAC_TIMESTEP Processes a single timestep for one CSAC.

    %% Unpack gating config
    absolute_flow_floor_kg_h = config.project.cutoff.flow_cutoff;
    delta_T_gate_threshold = config.project.cutoff.delta_T_gate_threshold;
    alpha_min = config.project.cutoff.alpha_min;

    %% ================================================================
    %% HEATING SEASON GATE
    %% ================================================================
    current_date = dateshift(time, 'start', 'day');
    [cs.season_state, actions] = manage_heating_season(...
        cs.season_state, current_date, t, daily_T_air_max_table, config, csac_id);

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
        return;
    end

    %% ================================================================
    %% MAIN PIPE TEMPERATURE ESTIMATION
    %% ================================================================
    num_active_houses = sum(current_data.flow_kg_h >= absolute_flow_floor_kg_h);
    is_csac_active = (num_active_houses >= config.project.initialization.min_active_houses);

    [T_junction_ukf_C, ukf_master_offset] = calculate_main_pipe_temp(...
        current_data, current_T_soil_C, cs.U_csac, absolute_flow_floor_kg_h, cs.ukf_states);

    if all(isnan(T_junction_ukf_C))
        return;
    end

    current_data.T_main_ukf_C = T_junction_ukf_C;
    cs.logger.timestamps(t) = time;

    %% ================================================================
    %% PER-HOUSE UKF UPDATES
    %% ================================================================
    gate_params_base = struct( ...
        'max_delta_T_change',     config.project.initialization.max_delta_T_change_rate, ...
        'hard_innovation_gate_C', config.project.initialization.max_innovation_C, ...
        'max_nis',                config.project.initialization.max_nis, ...
        'disable_stability_gate', config.project.debug.disable_stability_gate, ...
        'debug_print',            config.project.debug.print_update_result);

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
            gate_params = gate_params_base;
            gate_params.last_valid_T_main_C = cs.last_valid_T_main_C(i);
            gate_params.house_id = house_id;
            gate_params.time = time;

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

            % --- Consecutive rejection check ---
            [cs.ukf_states{i}, cs.consecutive_rejection_counters(i)] = ...
                check_consecutive_rejections(cs.ukf_states{i}, ...
                cs.consecutive_rejection_counters(i), ...
                update_result.rejected, update_result.accepted, config);
        end

        % --- Logging ---
        if log_ukf
            cs.logger = update_logger(cs.logger, t, i, time, cs.ukf_states{i}, update_result.diagnostics);
        elseif t > 1
            cs.logger.state_estimates(:, i, t) = cs.logger.state_estimates(:, i, t-1);
            cs.logger.covariance_posterior(:, i, t) = cs.logger.covariance_posterior(:, i, t-1);
        end
    end

    %% ================================================================
    %% SNAPSHOT UPDATE
    %% ================================================================
    [cs.ukf_states_snapshot, cs.season_state.snapshot_timestep] = ...
        update_snapshot(cs.season_state, current_date, t, ...
        cs.ukf_states, cs.ukf_states_snapshot, config);

    %% ================================================================
    %% MASTER OFFSET APPLICATION
    %% ================================================================
    cs.ukf_states = apply_master_offset(cs.ukf_states, ukf_master_offset, ...
        num_active_houses, config);

end