function [ukf_state, result] = update_house_ukf_gated(ukf_state, house_data, current_T_soil_C, config, gate_params)
%UPDATE_HOUSE_UKF_GATED Performs a gated UKF update for a single house.
%
%   This function wraps the UKF update with all pre- and post-update gates:
%   - Stability gate (checks that T_main hasn't changed too fast)
%   - Hard innovation gate (rejects large innovations before update)
%   - NIS gate (reverts update if normalized innovation squared is too high)
%   - Step size limiting is handled inside update_ukf_house itself
%
%   Args:
%       ukf_state (struct): Current UKF state for this house (.x, .P, .R, .Q).
%       house_data (table row): Single row of meter data for this house.
%           Must contain: T_supply_C, T_main_ukf_C, flow_kg_h, length_service_m.
%       current_T_soil_C (scalar): Current soil temperature in °C.
%       config (struct): Project configuration struct.
%       gate_params (struct): Gating parameters:
%           .last_valid_T_main_C (scalar): Previous valid T_main for stability check.
%           .max_delta_T_change (scalar): Max allowed T_main change rate.
%           .hard_innovation_gate_C (scalar): Max allowed absolute innovation.
%           .max_nis (scalar): Max allowed normalized innovation squared.
%           .disable_stability_gate (logical): If true, skip stability check.
%           .debug_print (logical): If true, print detailed update results.
%           .house_id (int): House ID for log messages.
%           .time (datetime): Current timestamp for log messages.
%
%   Returns:
%       ukf_state (struct): Updated (or reverted) UKF state.
%       result (struct): Outcome of the update attempt:
%           .accepted (logical): Whether the update was accepted.
%           .rejected (logical): Whether the update was rejected by a gate.
%           .skipped (logical): Whether the update was skipped (stability gate).
%           .diagnostics (struct): Diagnostics from update_ukf_house (if run).
%           .innovation_gate (scalar): Updated innovation gate value (if accepted).
%           .reason (char): Description of outcome.

    %% Initialize result
    result.accepted = false;
    result.rejected = false;
    result.skipped = false;
    result.diagnostics = struct();
    result.innovation_gate = NaN;
    result.reason = '';

    house_id = gate_params.house_id;
    time = gate_params.time;

    %% 1. Stability Gate
    if isnan(gate_params.last_valid_T_main_C)
        delta_T_main_change = 0;
    else
        delta_T_main_change = abs(house_data.T_main_ukf_C - gate_params.last_valid_T_main_C);
    end

    if gate_params.disable_stability_gate
        is_system_stable = true;
    else
        is_system_stable = (delta_T_main_change < gate_params.max_delta_T_change);
    end

    if ~is_system_stable
        result.skipped = true;
        result.reason = 'stability_gate';
        return;
    end

    %% 2. Pre-Update Innovation Gate
    predicted_temp = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, ...
        ukf_state.x(2), house_data.length_service_m, current_T_soil_C);
    innovation = house_data.T_supply_C - (predicted_temp - ukf_state.x(1));

    if abs(innovation) > gate_params.hard_innovation_gate_C
        result.rejected = true;
        result.reason = 'hard_innovation_gate';
        if gate_params.debug_print
            fprintf('UKF GATE REJECT | house %d | time %s | innovation=%.3f°C > gate=%.1f°C\n', ...
                house_id, string(time), innovation, gate_params.hard_innovation_gate_C);
        end
        return;
    end

    %% 3. Run UKF Update
    x_before = ukf_state.x;
    P_before = ukf_state.P;

    [ukf_state, diagnostics] = update_ukf_house(ukf_state, house_data, current_T_soil_C, config);

    % Check for NaN output
    if any(isnan(ukf_state.x), 'all') || any(isnan(ukf_state.P), 'all')
        warning('UKF returned NaN for house %d at %s — reverting', house_id, string(time));
        ukf_state.x = x_before;
        ukf_state.P = P_before;
        result.rejected = true;
        result.reason = 'nan_output';
        return;
    end

    %% 4. NIS Gate (post-update quality check)
    nis_value = NaN;
    if isfield(diagnostics, 'y') && isfield(diagnostics, 'P_zz') ...
            && isfinite(diagnostics.P_zz) && diagnostics.P_zz > 0
        nis_value = (diagnostics.y^2) / diagnostics.P_zz;
    end

    if isfinite(nis_value) && nis_value > gate_params.max_nis
        % Revert state
        ukf_state.x = x_before;
        ukf_state.P = P_before;
        result.rejected = true;
        result.reason = 'nis_gate';
        result.diagnostics = diagnostics;
        if gate_params.debug_print
            fprintf('UKF NIS REJECT | house %d | time %s | NIS=%.2f > %.1f | innov=%.3f\n', ...
                house_id, string(time), nis_value, gate_params.max_nis, diagnostics.y);
        end
        return;
    end

    %% 5. Update Accepted
    result.accepted = true;
    result.reason = 'accepted';
    result.diagnostics = diagnostics;

    % Update innovation gate from posterior
    if isfield(diagnostics, 'P_zz') && isfinite(diagnostics.P_zz) && diagnostics.P_zz > 0
        result.innovation_gate = max(1, sqrt(diagnostics.P_zz) * config.project.initialization.innovation_gate_N_sigma);
    end

    if gate_params.debug_print
        fprintf('UKF UPDATE | house %d | time %s | NIS=%.2f | innov=%.3f | dx=[%.4f %.6f]\n', ...
            house_id, string(time), nis_value, diagnostics.y, ...
            ukf_state.x(1) - x_before(1), ukf_state.x(2) - x_before(2));
    end
end