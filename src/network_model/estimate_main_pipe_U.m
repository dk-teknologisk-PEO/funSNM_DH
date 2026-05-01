function [U_main_new, diagnostics] = estimate_main_pipe_U(all_cs, csac_ids, topology, ...
    U_main_current, current_T_soil_C, config)
%ESTIMATE_MAIN_PIPE_U Estimates main pipe U-value by comparing fit quality
%   at slightly different U values.
%
%   Fits the network inlet temperature for U_current, U_current+dU, and
%   U_current-dU, then moves in the direction that reduces total residuals.

    cfg = config.project.main_pipe_U_estimation;
    num_csacs = numel(csac_ids);

    %% Extract CSAC inlet temperatures, positions, and flows
    positions = zeros(num_csacs, 1);
    T_inlet_measured = nan(num_csacs, 1);
    csac_flows = zeros(num_csacs, 1);
    csac_active = false(num_csacs, 1);

    for c = 1:num_csacs
        csac_topo_idx = find([topology.cul_de_sacs.id] == csac_ids(c));
        if ~isempty(csac_topo_idx)
            positions(c) = topology.cul_de_sacs(csac_topo_idx).dist_on_main_m;
        end

        cs = all_cs{c};
        if isfield(cs, 'T_inlet_fitted') && isfinite(cs.T_inlet_fitted) && cs.season_state.active
            T_inlet_measured(c) = cs.T_inlet_fitted;
            csac_active(c) = true;
        end
        if isfield(cs, 'current_total_flow')
            csac_flows(c) = cs.current_total_flow;
        end
    end

    %% Check minimum active CSACs
    if sum(csac_active) < 3
        U_main_new = U_main_current;
        diagnostics.adjusted = false;
        diagnostics.reason = 'insufficient_csacs';
        return;
    end

    %% Sort by position
    [positions_sorted, sort_idx] = sort(positions);
    T_measured_sorted = T_inlet_measured(sort_idx);
    flows_sorted = csac_flows(sort_idx);
    active_sorted = csac_active(sort_idx);

    Q_total = cumsum(flows_sorted, 'reverse');
    x_diff = positions_sorted;
    x_diff(2:end) = diff(positions_sorted);

    valid_idx = find(active_sorted);
    if numel(valid_idx) < 3
        U_main_new = U_main_current;
        diagnostics.adjusted = false;
        diagnostics.reason = 'insufficient_valid';
        return;
    end

    %% Evaluate fit quality at three U values
    dU = 0.005;
    U_test = [U_main_current - dU, U_main_current, U_main_current + dU];
    costs = zeros(1, 3);

    for u = 1:3
        costs(u) = fit_main_pipe_cost(U_test(u), T_measured_sorted, Q_total, ...
            x_diff, current_T_soil_C, valid_idx);
    end

    %% Determine direction of improvement
    % costs(1) = cost at U-dU, costs(2) = cost at U, costs(3) = cost at U+dU
    if costs(1) < costs(2) && costs(1) < costs(3)
        % Lower U is better
        direction = -1;
    elseif costs(3) < costs(2) && costs(3) < costs(1)
        % Higher U is better
        direction = 1;
    else
        % Current U is best (local minimum)
        U_main_new = U_main_current;
        diagnostics.adjusted = false;
        diagnostics.reason = 'at_minimum';
        diagnostics.costs = costs;
        diagnostics.U_test = U_test;
        return;
    end

    %% Compute adjustment magnitude from gradient
    % Numerical gradient of cost with respect to U
    cost_gradient = (costs(3) - costs(1)) / (2 * dU);

    % Step size proportional to gradient, but clamped
    raw_adjustment = -cfg.gain * cost_gradient;
    adjustment = max(-cfg.max_adjustment_per_step, min(cfg.max_adjustment_per_step, raw_adjustment));

    % Ensure adjustment is in the correct direction
    if sign(adjustment) ~= direction
        adjustment = direction * abs(adjustment);
    end

    U_main_new = U_main_current + adjustment;
    U_main_new = max(cfg.U_min, min(cfg.U_max, U_main_new));

    %% Diagnostics
    diagnostics.adjusted = true;
    diagnostics.reason = 'updated';
    diagnostics.costs = costs;
    diagnostics.U_test = U_test;
    diagnostics.cost_gradient = cost_gradient;
    diagnostics.direction = direction;
    diagnostics.adjustment = adjustment;
    diagnostics.U_old = U_main_current;
    diagnostics.U_new = U_main_new;
end


function cost = fit_main_pipe_cost(U_main, T_measured_sorted, Q_total, x_diff, T_soil, valid_idx)
%FIT_MAIN_PIPE_COST Fits T_network_inlet for a given U_main and returns
%   the weighted sum of squared residuals.

    % Optimize T_network_inlet for this U_main
    T_start_guess = max(T_measured_sorted(valid_idx)) + 1.0;

    cost_fun = @(T_start) compute_residual_cost(T_start, T_measured_sorted, ...
        Q_total, x_diff, T_soil, U_main, valid_idx);

    opts = optimset('Display', 'off', 'TolX', 0.01, 'TolFun', 0.01, 'MaxIter', 100);
    [~, cost] = fminsearch(cost_fun, T_start_guess, opts);
end


function cost = compute_residual_cost(T_start, T_measured, Q_total, x_diff, T_soil, U_main, valid_idx)
%COMPUTE_RESIDUAL_COST Forward-simulates main pipe and computes residuals.
    num = numel(T_measured);
    T_sim = zeros(num, 1);
    T_prev = T_start;
    for c = 1:num
        T_prev = get_supply_temp(T_prev, Q_total(c), U_main, x_diff(c), T_soil);
        T_sim(c) = T_prev;
    end
    residuals = T_sim(valid_idx) - T_measured(valid_idx);
    cost = sum(residuals.^2);
end