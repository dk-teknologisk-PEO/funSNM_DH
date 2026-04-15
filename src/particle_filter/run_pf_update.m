function csac_state = run_pf_update(csac_state, i, house_data, current_T_soil_C, params)
%RUN_PF_UPDATE Orchestrate PF update for one house with all gates.
%   Applies stability gate, innovation gate, NIS gate, then calls
%   update_pf_house for the core PF math. Updates csac_state in place.
%
%   This is the gating/orchestration layer. The core PF math lives in
%   update_pf_house.m.

    house_id = csac_state.house_ids(i);
    
    %% Stability gate
    if isnan(csac_state.last_valid_T_main_pf_C(i))
        delta_T_main_change = 0;
    else
        delta_T_main_change = abs(house_data.T_main_pf_C - csac_state.last_valid_T_main_pf_C(i));
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
    predicted_temp = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
        csac_state.pf_states{i}.x(2), house_data.length_service_m, current_T_soil_C);
    innovation = house_data.T_supply_C - (predicted_temp - csac_state.pf_states{i}.x(1));
    
    if abs(innovation) > params.hard_innovation_gate_C
        csac_state.gate_reject_count_pf = csac_state.gate_reject_count_pf + 1;
        if params.debug_print_update_result
            fprintf('PF  GATE REJECT | house %d | innovation=%.3f°C > gate=%.1f°C\n', ...
                house_id, innovation, params.hard_innovation_gate_C);
        end
        return;
    end
    
    %% Core PF update
    x_before = csac_state.pf_states{i}.x;
    
    [csac_state.pf_particles{i}, est_pf, cov_pf, diagnostics_pf] = update_pf_house(...
        csac_state.pf_particles{i}, house_data, current_T_soil_C, ...
        params.R_base, params.Q_base, params.config);
    
    if any(isnan(est_pf), 'all') || any(isnan(cov_pf), 'all')
        warning('PF returned NaN for house %d', house_id);
        return;
    end
    
    %% NIS gate (post-update check)
    nis_value = NaN;
    if isfield(diagnostics_pf, 'P_zz') && isfinite(diagnostics_pf.P_zz) && diagnostics_pf.P_zz > 0
        if isfield(diagnostics_pf, 'y') && isfinite(diagnostics_pf.y)
            nis_value = (diagnostics_pf.y^2) / diagnostics_pf.P_zz;
        end
    end
    
    if isfinite(nis_value) && nis_value > params.max_nis
        csac_state.pf_states{i}.x = x_before;
        csac_state.gate_reject_count_pf = csac_state.gate_reject_count_pf + 1;
        if params.debug_print_update_result
            fprintf('PF  NIS REJECT | house %d | NIS=%.2f > %.1f | innov=%.3f\n', ...
                house_id, nis_value, params.max_nis, diagnostics_pf.y);
        end
        return;
    end
    
    %% Accept update — record results
    csac_state.gate_accept_count_pf = csac_state.gate_accept_count_pf + 1;
    csac_state.pf_states{i}.x = est_pf;
    csac_state.pf_states{i}.P = cov_pf;
    csac_state.pf_offsets(i) = csac_state.pf_states{i}.x(1);
    csac_state.last_valid_T_main_pf_C(i) = house_data.T_main_pf_C;
    csac_state.pf_updated(i) = true;
    csac_state.diagnostics_pf{i} = diagnostics_pf;
    
    if params.debug_print_update_result
        fprintf('PF  UPDATE | house %d | NIS=%.2f | innov=%.3f | dx=[%.4f %.6f]\n', ...
            house_id, nis_value, diagnostics_pf.y, ...
            csac_state.pf_states{i}.x(1) - x_before(1), ...
            csac_state.pf_states{i}.x(2) - x_before(2));
    end
end