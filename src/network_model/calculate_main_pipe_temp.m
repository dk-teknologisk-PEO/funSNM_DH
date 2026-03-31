function T_main_fit_C = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_cutoff, states, T_start_C)
% CALCULATE_MAIN_PIPE_TEMP Estimates the main pipe temperature at each
% house junction using a leave-one-out optimization approach.
%
% IMPORTANT: 'states' must be aligned with 'current_data' rows, i.e.,
% states{i} corresponds to current_data(i,:). The caller is responsible
% for ensuring this alignment.
    
    num_houses = size(current_data,1);
    T_main_C = nan([num_houses,1]);

    T_main_fit_C = nan([num_houses,1]);

    % segment lengths along main pipe
    x_diff_m = [current_data.x_pos_m];
    x_diff_m(2:end) = diff(x_diff_m);
    
    % cumulative flow (from last house backward = total flow at each
    % segment)
    Q_total = cumsum(current_data.flow_kg_h,'reverse');
    
    if nargin == 6 && ~isempty(T_start_C) && ~isnan(T_start_C)
        % Direct forward simulation from known inlet temperature
        T_prev_C = T_start_C;
        for i=1:num_houses
            T_prev_C = get_supply_temp(T_prev_C, Q_total(i), U_csac, x_diff_m(i),current_T_soil_C);
            T_main_C(i)=T_prev_C;
        end
        return
    end

    % --- Estimate T_main at each junction from state estimates ---
    % Extract current state estimates (must be aligned with current_data rows)
    u_service = cellfun(@(s) s.x(2), states(:));  % force column
    offsets   = cellfun(@(s) s.x(1), states(:));

    % sanity check: dimensions must match
    if length(u_service) ~= num_houses
        warning('calculate_main_pipe_temp: states (%d) and current_data (%d) size mismatch.', ...
            length(u_service), num_houses);
        return
    end

    service_lengths_m = current_data.length_service_m;

    % Back-calculate what T_main must have been, given measured T_supply and state estimates
    T_main_C = get_main_temp( ...
        current_data.T_supply_C + offsets, ...
        current_data.flow_kg_h, ...
        u_service, ...
        service_lengths_m, ...
        current_T_soil_C);

    % Protect against Inf/NaN from low-flow divisions
    T_main_C(~isfinite(T_main_C)) = NaN;
    current_data.T_main_C = T_main_C;

    % Robust starting point for fminsearch
    valid_T_main = T_main_C(isfinite(T_main_C));
    if isempty(valid_T_main)
        return  % Cannot estimate anything
    end
    default_start = max(valid_T_main);

    % --- Leave-one-out optimization for each house ---
    opts = optimset('Display', 'off', 'TolX', 0.05, 'TolFun', 0.05, 'MaxIter', 200);

    for m = 1:num_houses
        if current_data.flow_kg_h(m) < flow_cutoff
            continue  % T_main_fit_C(m) stays NaN
        end

        current_data_sub = current_data;
        current_data_sub(m, :) = [];

        % Check that we have enough valid reference houses
        n_valid_ref = sum(isfinite(current_data_sub.T_main_C));
        if n_valid_ref < 2
            continue  % Not enough data for a meaningful fit
        end

        fun = @(x) main_temp_estimator(x, current_data, current_data_sub, current_T_soil_C, U_csac, NaN);

        start_val = default_start;
        if ~isfinite(start_val)
            continue
        end

        [p, ~, exitflag] = fminsearch(fun, start_val, opts);

        if exitflag <= 0 || ~isfinite(p)
            continue
        end

        % Forward-simulate from optimized inlet to get T_main at house m
        T_prev_C = p(1);
        for i = 1:m
            T_prev_C = get_supply_temp(T_prev_C, Q_total(i), U_csac, x_diff_m(i), current_T_soil_C);
        end
        T_main_fit_C(m) = T_prev_C;
    end
end