function csac_state = process_timestep(csac_state, t, time, T_soil_C, T_air_C, ...
        daily_T_air_max_table, params)
%PROCESS_TIMESTEP Process one timestep for a CSAC.

    current_data = csac_state.meter_data_csac(csac_state.meter_data_csac.timestamp == time, :);
    current_data = sortrows(current_data, 'house_id');
    
    current_T_soil_C = T_soil_C(T_soil_C.time == time, :).values;
    
    if isempty(current_data) || isempty(current_T_soil_C)
        return;
    end
    
    %% Season gate
    current_date = dateshift(time, 'start', 'day');
    csac_state = evaluate_season_transition(csac_state, t, current_date, ...
        daily_T_air_max_table, params);
    
    %% If not active, carry forward and return
    if ~csac_state.season_active
        if t > 1
            for i = 1:csac_state.num_houses
                csac_state.logger_ukf.state_estimates(:, i, t) = csac_state.logger_ukf.state_estimates(:, i, t-1);
                csac_state.logger_ukf.covariance_posterior(:, i, t) = csac_state.logger_ukf.covariance_posterior(:, i, t-1);
                csac_state.logger_pf.state_estimates(:, i, t) = csac_state.logger_pf.state_estimates(:, i, t-1);
                csac_state.logger_pf.covariance_posterior(:, i, t) = csac_state.logger_pf.covariance_posterior(:, i, t-1);
            end
        end
        csac_state.logger_ukf.timestamps(t) = time;
        csac_state.logger_pf.timestamps(t) = time;
        return;
    end
    
    %% Active season processing
    num_active_houses = sum(current_data.flow_kg_h >= params.absolute_flow_floor_kg_h);
    is_csac_active = (num_active_houses >= params.min_active_houses);
    
    [T_junction_ukf_C, ukf_master_offset] = calculate_main_pipe_temp(...
        current_data, current_T_soil_C, csac_state.U_csac, ...
        params.absolute_flow_floor_kg_h, csac_state.ukf_states);
    [T_junction_pf_C, pf_master_offset] = calculate_main_pipe_temp(...
        current_data, current_T_soil_C, csac_state.U_csac, ...
        params.absolute_flow_floor_kg_h, csac_state.pf_states);
    
    if all(isnan(T_junction_ukf_C)) || all(isnan(T_junction_pf_C))
        return;
    end
    
    current_data.T_main_ukf_C = T_junction_ukf_C;
    current_data.T_main_pf_C = T_junction_pf_C;
    csac_state.logger_ukf.timestamps(t) = time;
    csac_state.logger_pf.timestamps(t) = time;
    
    % Reset per-timestep flags
    csac_state.ukf_updated = false(1, csac_state.num_houses);
    csac_state.diagnostics_ukf = cell(1, csac_state.num_houses);
    csac_state.pf_updated = false(1, csac_state.num_houses);
    csac_state.diagnostics_pf = cell(1, csac_state.num_houses);
    
    %% Per-house updates
    for i = 1:csac_state.num_houses
        house_data = current_data(current_data.house_id == csac_state.house_ids(i), :);
        if isempty(house_data), continue; end
        
        % Flow and delta-T gates
        [is_flow_ok_ukf, ~, ~, ~] = flow_validity_gate(house_data.flow_kg_h, ...
            csac_state.ukf_states{i}.x(2), house_data.length_service_m, params.alpha_min);
        is_flow_ok_ukf = is_flow_ok_ukf && (house_data.flow_kg_h > params.absolute_flow_floor_kg_h);
        delta_T_ok_ukf = (house_data.T_main_ukf_C - current_T_soil_C) >= params.delta_T_gate_threshold;
        
        [is_flow_ok_pf, ~, ~, ~] = flow_validity_gate(house_data.flow_kg_h, ...
            csac_state.pf_states{i}.x(2), house_data.length_service_m, params.alpha_min);
        is_flow_ok_pf = is_flow_ok_pf && (house_data.flow_kg_h > params.absolute_flow_floor_kg_h);
        delta_T_ok_pf = (house_data.T_main_pf_C - current_T_soil_C) >= params.delta_T_gate_threshold;
        
        % UKF update
        if is_csac_active && is_flow_ok_ukf && delta_T_ok_ukf
            csac_state = run_ukf_update(csac_state, i, house_data, current_T_soil_C, params);
        end
        
        % PF update
        if is_csac_active && is_flow_ok_pf && delta_T_ok_pf
            csac_state = run_pf_update(csac_state, i, house_data, current_T_soil_C, params);
        end
        
        % Logger update
        if csac_state.ukf_updated(i)
            csac_state.logger_ukf = update_logger(csac_state.logger_ukf, t, i, time, ...
                csac_state.ukf_states{i}, csac_state.diagnostics_ukf{i});
        elseif t > 1
            csac_state.logger_ukf.state_estimates(:, i, t) = csac_state.logger_ukf.state_estimates(:, i, t-1);
            csac_state.logger_ukf.covariance_posterior(:, i, t) = csac_state.logger_ukf.covariance_posterior(:, i, t-1);
        end
        % Same for PF (using csac_state.pf_updated flag)
        if csac_state.pf_updated(i)
            csac_state.logger_pf = update_logger(csac_state.logger_pf, t, i, time, ...
                csac_state.pf_states{i}, csac_state.diagnostics_pf{i});
        elseif t > 1
            csac_state.logger_pf.state_estimates(:, i, t) = csac_state.logger_pf.state_estimates(:, i, t-1);
            csac_state.logger_pf.covariance_posterior(:, i, t) = csac_state.logger_pf.covariance_posterior(:, i, t-1);
        end
    end
    
    %% Snapshot update
    csac_state.snapshot_timestep = t;
    for i = 1:csac_state.num_houses
        csac_state.ukf_states_snapshot{i}.x = csac_state.ukf_states{i}.x;
        csac_state.ukf_states_snapshot{i}.P = csac_state.ukf_states{i}.P;
        csac_state.pf_states_snapshot{i}.x = csac_state.pf_states{i}.x;
        csac_state.pf_states_snapshot{i}.P = csac_state.pf_states{i}.P;
    end
    
    %% Master offset
    csac_state = apply_master_offset(csac_state, ukf_master_offset, pf_master_offset, ...
        num_active_houses, params);
end