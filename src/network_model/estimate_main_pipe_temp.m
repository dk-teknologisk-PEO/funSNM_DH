function [T_csac_inlet_C, main_pipe_diag] = estimate_main_pipe_temp(...
    all_cs, csac_ids, topology, current_T_soil_C, U_main, config)
%ESTIMATE_MAIN_PIPE_TEMP Estimates main pipe temperature profile from CSAC data.
%
%   Uses each CSAC's fitted inlet temperature as a "measurement" of the
%   main pipe temperature at that CSAC's junction. Fits the network inlet
%   temperature by minimizing weighted residuals.
%
%   This is the main-pipe equivalent of calculate_main_pipe_temp.
%
%   Args:
%       all_cs (cell array): CSAC state structs (must have .T_inlet_fitted set).
%       csac_ids (int array): CSAC identifiers.
%       topology (struct): Network topology with cul_de_sacs positions.
%       current_T_soil_C (scalar): Current soil temperature.
%       U_main (scalar): Main pipe U-value [W/m/K].
%       config (struct): Project configuration.
%
%   Returns:
%       T_csac_inlet_C (Nx1 double): Estimated temperature at each CSAC junction.
%       main_pipe_diag (struct): Diagnostic information.

    num_csacs = numel(csac_ids);
    T_csac_inlet_C = nan(num_csacs, 1);

    %% Extract CSAC junction positions and flows
    csac_positions = zeros(num_csacs, 1);
    csac_T_inlet_measured = nan(num_csacs, 1);
    csac_total_flow = zeros(num_csacs, 1);
    csac_weights = zeros(num_csacs, 1);
    csac_active = false(num_csacs, 1);

    for c = 1:num_csacs
        % Position on main pipe
        csac_topo_idx = find([topology.cul_de_sacs.id] == csac_ids(c));
        if ~isempty(csac_topo_idx)
            csac_positions(c) = topology.cul_de_sacs(csac_topo_idx).dist_on_main_m;
        end

        % Check if this CSAC has a valid T_inlet estimate
        cs = all_cs{c};
        if isfield(cs, 'T_inlet_fitted') && isfinite(cs.T_inlet_fitted)
            csac_T_inlet_measured(c) = cs.T_inlet_fitted;
            csac_active(c) = true;
        end

        % Total flow for this CSAC (sum of all house flows at current timestep)
        if isfield(cs, 'current_total_flow')
            csac_total_flow(c) = cs.current_total_flow;
        end
    end

    %% Check if we have enough active CSACs
    if sum(csac_active) < 2
        main_pipe_diag.valid = false;
        main_pipe_diag.reason = 'insufficient_csacs';
        return;
    end

    %% Sort by position (should already be sorted, but ensure it)
    [csac_positions_sorted, sort_idx] = sort(csac_positions);
    csac_T_measured_sorted = csac_T_inlet_measured(sort_idx);
    csac_flow_sorted = csac_total_flow(sort_idx);
    csac_active_sorted = csac_active(sort_idx);

    %% Compute cumulative flow (sum of all downstream CSAC flows)
    Q_total = cumsum(csac_flow_sorted, 'reverse');

    %% Compute segment distances
    x_diff_m = csac_positions_sorted;
    x_diff_m(2:end) = diff(csac_positions_sorted);

    %% Compute weights for active CSACs
    % Weight by number of active houses (more houses = more reliable T_inlet)
    for c = 1:num_csacs
        orig_idx = sort_idx(c);
        if csac_active_sorted(c)
            cs = all_cs{orig_idx};
            csac_weights(c) = cs.gate_accept_count + 1;  % +1 to avoid zero
        end
    end

    valid_csacs = find(csac_active_sorted & csac_weights > 0);
    if numel(valid_csacs) < 2
        main_pipe_diag.valid = false;
        main_pipe_diag.reason = 'insufficient_valid_csacs';
        return;
    end

    fit_weights = csac_weights(valid_csacs);
    fit_weights = fit_weights / sum(fit_weights);

    %% Optimize: find T_network_inlet that best fits the CSAC inlet temperatures
    opts = optimset('Display', 'off', 'TolX', 0.05, 'TolFun', 0.05, 'MaxIter', 200);

    % Initial guess: highest measured T_inlet (closest to source)
    T_start_guess = max(csac_T_measured_sorted(csac_active_sorted));

    cost_fun = @(T_start) main_pipe_cost(T_start, csac_T_measured_sorted, ...
        Q_total, x_diff_m, current_T_soil_C, U_main, valid_csacs, fit_weights);

    [T_network_inlet, ~, exitflag] = fminsearch(cost_fun, T_start_guess, opts);

    if exitflag <= 0 || ~isfinite(T_network_inlet)
        main_pipe_diag.valid = false;
        main_pipe_diag.reason = 'optimization_failed';
        return;
    end

    %% Forward-simulate main pipe temperatures from fitted inlet
    T_prev = T_network_inlet;
    T_csac_fitted_sorted = zeros(num_csacs, 1);
    for c = 1:num_csacs
        T_prev = get_supply_temp(T_prev, Q_total(c), U_main, x_diff_m(c), current_T_soil_C);
        T_csac_fitted_sorted(c) = T_prev;
    end

    %% Unsort back to original CSAC order
    T_csac_inlet_C(sort_idx) = T_csac_fitted_sorted;

    %% Diagnostics
    main_pipe_diag.valid = true;
    main_pipe_diag.T_network_inlet = T_network_inlet;
    main_pipe_diag.T_csac_fitted = T_csac_inlet_C;
    main_pipe_diag.T_csac_measured = csac_T_inlet_measured;
    main_pipe_diag.residuals = T_csac_inlet_C - csac_T_inlet_measured;
    main_pipe_diag.csac_positions = csac_positions;
end


function cost = main_pipe_cost(T_start, T_measured, Q_total, x_diff, T_soil, U_main, valid_idx, weights)
%MAIN_PIPE_COST Cost function for main pipe inlet temperature optimization.
    num = numel(T_measured);
    T_sim = zeros(num, 1);
    T_prev = T_start;
    for c = 1:num
        T_prev = get_supply_temp(T_prev, Q_total(c), U_main, x_diff(c), T_soil);
        T_sim(c) = T_prev;
    end
    residuals = T_sim(valid_idx) - T_measured(valid_idx);
    cost = sum(weights .* residuals.^2);
end