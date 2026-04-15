function csac_state = run_ukf_update(csac_state, i, house_data, current_T_soil_C, params)
%RUN_UKF_UPDATE Orchestrate UKF update for one house with all gates.
%   Applies stability gate, innovation gate, NIS gate, then calls
%   update_ukf_house for the core UKF math. Updates csac_state in place.
%
%   This is the gating/orchestration layer. The core UKF math lives in
%   update_ukf_house.m.

    house_id = csac_state.house_ids(i);
    
    %% Stability gate
    if isnan(csac_state.last_valid_T_main_ukf_C(i))
        delta_T_main_change = 0;
    else
        delta_T_main_change = abs(house_data.T_main_ukf_C - csac_state.last_valid_T_main_ukf_C(i));
    end
    
    if params.debug_disable_stability_gate
        is_system_stable = true;
    else
        is_system_stable = (delta_T_main_change < params.max_delta_T_change);
    end
    
    if ~is_system_stable
        return;
    end
    
    %% Innovation gate (pre-update check)
    predicted_temp = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, ...
        csac_state.ukf_states{i}.x(2), house_data.length_service_m, current_T_soil_C);
    innovation = house_data.T_supply_C - (predicted_temp - csac_state.ukf_states{i}.x(1));
    
    if abs(innovation) > params.hard_innovation_gate_C
        csac_state.gate_reject_count_ukf = csac_state.gate_reject_count_ukf + 1;
        if params.debug_print_update_result
            fprintf('UKF GATE REJECT | house %d | innovation=%.3f°C > gate=%.1f°C\n', ...
                house_id, innovation, params.hard_innovation_gate_C);
        end
        return;
    end
    
    %% Core UKF update
    x_before = csac_state.ukf_states{i}.x;
    P_before = csac_state.ukf_states{i}.P;
    
    [csac_state.ukf_states{i}, diagnostics_ukf] = update_ukf_house(...
        csac_state.ukf_states{i}, house_data, current_T_soil_C, params.config);
    
    if any(isnan(csac_state.ukf_states{i}.x), 'all') || ...
            any(isnan(csac_state.ukf_states{i}.P), 'all')
        warning('UKF returned NaN for house %d', house_id);
        csac_state.ukf_states{i}.x = x_before;
        csac_state.ukf_states{i}.P = P_before;
        return;
    end
    
    %% NIS gate (post-update check)
    nis_value = NaN;
    if isfield(diagnostics_ukf, 'y') && isfield(diagnostics_ukf, 'P_zz') ...
            && isfinite(diagnostics_ukf.P_zz) && diagnostics_ukf.P_zz > 0
        nis_value = (diagnostics_ukf.y^2) / diagnostics_ukf.P_zz;
    end
    
    if isfinite(nis_value) && nis_value > params.max_nis
        csac_state.ukf_states{i}.x = x_before;
        csac_state.ukf_states{i}.P = P_before;
        csac_state.gate_reject_count_ukf = csac_state.gate_reject_count_ukf + 1;
        if params.debug_print_update_result
            fprintf('UKF NIS REJECT | house %d | NIS=%.2f > %.1f | innov=%.3f\n', ...
                house_id, nis_value, params.max_nis, diagnostics_ukf.y);
        end
        return;
    end
    
    %% Accept update — record results
    csac_state.gate_accept_count_ukf = csac_state.gate_accept_count_ukf + 1;
    csac_state.ukf_offsets(i) = csac_state.ukf_states{i}.x(1);
    csac_state.last_valid_T_main_ukf_C(i) = house_data.T_main_ukf_C;
    csac_state.last_update_timestamp_ukf(i) = house_data.timestamp;
    csac_state.ukf_updated(i) = true;
    csac_state.diagnostics_ukf{i} = diagnostics_ukf;
    
    if params.debug_print_update_result
        fprintf('UKF UPDATE | house %d | NIS=%.2f | innov=%.3f | dx=[%.4f %.6f]\n', ...
            house_id, nis_value, diagnostics_ukf.y, ...
            csac_state.ukf_states{i}.x(1) - x_before(1), ...
            csac_state.ukf_states{i}.x(2) - x_before(2));
    end
end