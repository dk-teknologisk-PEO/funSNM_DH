function [T_main_fit_C, master_offset_C] = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_cutoff, states, T_start_C)
% CALCULATE_MAIN_PIPE_TEMP Estimates main pipe temperature and a master offset.
% A single optimization is performed to find the best inlet temperature and
% a shared meter offset for the entire cul-de-sac.
%
% Houses are weighted by flow quality — high-flow houses contribute more
% to the fit than low-flow houses near the cutoff threshold.
%
% If T_start_C is provided (measured inlet temperature), the optimization
% is skipped and T_main is forward-simulated directly from that value.

    num_houses = size(current_data,1);
    
    % Initialize outputs
    T_main_fit_C = nan(num_houses,1);
    master_offset_C = NaN;
    
    % Pre-compute segment geometry
    Q_total = cumsum(current_data.flow_kg_h, 'reverse');
    x_diff_m = current_data.x_pos_m;
    x_diff_m(2:end) = diff(x_diff_m);
    
    % --- Check if we have a measured inlet temperature ---
    has_measured_T_start = (nargin >= 6) && ~isempty(T_start_C) && isfinite(T_start_C);
    
    if has_measured_T_start
        % === DIRECT PATH: forward-simulate from measured T_start ===
        master_offset_C = 0;
        T_prev_C = T_start_C;
        for i = 1:num_houses
            T_prev_C = get_supply_temp(T_prev_C, Q_total(i), U_csac, x_diff_m(i), current_T_soil_C);
            if current_data.flow_kg_h(i) >= flow_cutoff
                T_main_fit_C(i) = T_prev_C;
            end
        end
        return;
    end
    
    % === OPTIMIZATION PATH: estimate T_start and master_offset ===
    
    % Back-calculate T_main_C from raw meter data
    u_service = cellfun(@(s) s.x(2), states(:));
    service_lengths_m = current_data.length_service_m;
    T_main_C = get_main_temp(current_data.T_supply_C, current_data.flow_kg_h, u_service, service_lengths_m, current_T_soil_C);
    T_main_C(~isfinite(T_main_C)) = NaN;
    current_data.T_main_C = T_main_C;
    
    % --- Compute flow-based weights for each house ---
    % Higher flow = more reliable T_main back-calculation
    c = 4186;
    weights = zeros(num_houses, 1);
    for i = 1:num_houses
        if current_data.flow_kg_h(i) >= flow_cutoff && isfinite(T_main_C(i))
            flow_kg_s = current_data.flow_kg_h(i) / 3600;
            theta = (u_service(i) * service_lengths_m(i)) / (2 * c);
            if flow_kg_s > theta
                alpha = (flow_kg_s - theta) / (flow_kg_s + theta);
            else
                alpha = 0;
            end
            % Weight is alpha squared — strongly favors high-flow houses
            weights(i) = alpha^2;
        end
    end
    
    % Global optimization
    valid_houses_for_fit = find(weights > 0);
    if numel(valid_houses_for_fit) < 2
        return;
    end
    
    % Normalize weights
    fit_weights = weights(valid_houses_for_fit);

    % Normalize weights
    % Cap maximum weight to prevent any single house from dominating
    max_weight_fraction = 0.35;  % No house contributes more than 35%
    fit_weights = min(fit_weights, max_weight_fraction * sum(fit_weights));
    
    % Re-normalize after capping
    fit_weights = fit_weights / sum(fit_weights);
    
    opts = optimset('Display', 'off', 'TolX', 0.05, 'TolFun', 0.05, 'MaxIter', 200);
    
    fun = @(params) main_temp_estimator_global_weighted(params, current_data, valid_houses_for_fit, fit_weights, current_T_soil_C, U_csac);
    default_start_T = max(T_main_C(valid_houses_for_fit));
    start_vector = [default_start_T, 0];
    
    [p_global, ~, exitflag] = fminsearch(fun, start_vector, opts);
    
    if exitflag <= 0 || any(~isfinite(p_global))
        warning('Global T_main optimization failed.');
        return;
    end
    
    optimized_T_start = p_global(1);
    master_offset_C = p_global(2);
    
    % Forward-simulate for all houses from the optimized T_start
    T_prev_C = optimized_T_start;
    for i = 1:num_houses
        T_prev_C = get_supply_temp(T_prev_C, Q_total(i), U_csac, x_diff_m(i), current_T_soil_C);
        if current_data.flow_kg_h(i) >= flow_cutoff
            T_main_fit_C(i) = T_prev_C;
        end
    end
end

