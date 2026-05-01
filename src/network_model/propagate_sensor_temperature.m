function [T_junctions, T_uncertainties] = propagate_sensor_temperature(...
    sensor_T, sensor_pos, sensor_noise_std, ...
    junction_positions, junction_flows, ...
    U_main, U_main_uncertainty, T_soil, topology)
%PROPAGATE_SENSOR_TEMPERATURE Propagates a reference sensor temperature
%   to all CSAC junctions on the main pipe, with distance-dependent uncertainty.
%
%   Propagates both upstream (toward plant) and downstream (away from plant)
%   from the sensor location.
%
%   Args:
%       sensor_T (scalar): Measured temperature at sensor [°C].
%       sensor_pos (scalar): Position of sensor on main pipe [m].
%       sensor_noise_std (scalar): Sensor measurement noise std [°C].
%       junction_positions (Nx1 double): Position of each CSAC junction [m].
%       junction_flows (Nx1 double): Total flow at each CSAC junction [kg/h].
%       U_main (scalar): Current main pipe U-value estimate [W/m/K].
%       U_main_uncertainty (scalar): Uncertainty in U_main [W/m/K].
%       T_soil (scalar): Current soil temperature [°C].
%       topology (struct): Network topology.
%
%   Returns:
%       T_junctions (Nx1 double): Estimated temperature at each junction.
%       T_uncertainties (Nx1 double): Uncertainty (std) at each junction.

    num_junctions = numel(junction_positions);
    T_junctions = nan(num_junctions, 1);
    T_uncertainties = nan(num_junctions, 1);

    % Sort junctions by position
    [pos_sorted, sort_idx] = sort(junction_positions);
    flows_sorted = junction_flows(sort_idx);

    % Cumulative flow (sum of all downstream)
    Q_total = cumsum(flows_sorted, 'reverse');

    % Segment distances
    x_diff = pos_sorted;
    x_diff(2:end) = diff(pos_sorted);

    % Find which sorted index the sensor is closest to
    % The sensor is AT a junction (or at the plant)
    [~, sensor_sorted_idx] = min(abs(pos_sorted - sensor_pos));

    % Check if sensor is at a junction or at the plant (before all junctions)
    sensor_is_at_plant = (sensor_pos < pos_sorted(1));

    %% Forward propagation (downstream from sensor)
    if sensor_is_at_plant
        start_idx = 1;
        T_prev = sensor_T;
        sigma_T_prev = sensor_noise_std;
    else
        start_idx = sensor_sorted_idx + 1;
        T_prev = sensor_T;
        sigma_T_prev = sensor_noise_std;

        % Set the sensor junction
        T_junctions(sort_idx(sensor_sorted_idx)) = sensor_T;
        T_uncertainties(sort_idx(sensor_sorted_idx)) = sensor_noise_std;
    end

    % Forward: sensor → downstream junctions
    for c = start_idx:num_junctions
        if sensor_is_at_plant
            dist = pos_sorted(c);
            if c == 1
                seg_dist = pos_sorted(c);
            else
                seg_dist = x_diff(c);
            end
        else
            seg_dist = x_diff(c);
        end

        % Propagate temperature
        T_out = get_supply_temp(T_prev, Q_total(c), U_main, seg_dist, T_soil);

        % Propagate uncertainty
        % 1. Temperature uncertainty propagation through alpha
        c_water = 4186;
        theta = (U_main * seg_dist) / (2 * c_water);
        flow_kg_s = Q_total(c) / 3600;
        if flow_kg_s > theta
            alpha_val = (flow_kg_s - theta) / (flow_kg_s + theta);
        else
            alpha_val = 0;
        end
        sigma_T_from_input = abs(alpha_val) * sigma_T_prev;

        % 2. U-value uncertainty contribution
        dU = 0.001;
        T_out_perturbed = get_supply_temp(T_prev, Q_total(c), U_main + dU, seg_dist, T_soil);
        dTdU = abs(T_out_perturbed - T_out) / dU;
        sigma_T_from_U = dTdU * U_main_uncertainty;

        % Combined uncertainty
        sigma_T_out = sqrt(sigma_T_from_input^2 + sigma_T_from_U^2);

        T_junctions(sort_idx(c)) = T_out;
        T_uncertainties(sort_idx(c)) = sigma_T_out;

        T_prev = T_out;
        sigma_T_prev = sigma_T_out;
    end

    %% Backward propagation (upstream from sensor)
    if ~sensor_is_at_plant && sensor_sorted_idx > 1
        T_prev = sensor_T;
        sigma_T_prev = sensor_noise_std;

        for c = (sensor_sorted_idx - 1):-1:1
            % Distance from junction c to junction c+1
            seg_dist = pos_sorted(c + 1) - pos_sorted(c);

            % Back-calculate: T_in from T_out
            % T_out = alpha * T_in + (1-alpha) * T_soil
            % T_in = (T_out - (1-alpha) * T_soil) / alpha
            c_water = 4186;
            theta = (U_main * seg_dist) / (2 * c_water);
            flow_kg_s = Q_total(c + 1) / 3600;
            if flow_kg_s > theta
                alpha_val = (flow_kg_s - theta) / (flow_kg_s + theta);
            else
                alpha_val = 0.01; % avoid division by zero
            end

            T_in = (T_prev - (1 - alpha_val) * T_soil) / alpha_val;

            % Uncertainty propagation (backward)
            sigma_T_from_output = sigma_T_prev / abs(alpha_val);

            % U-value uncertainty
            dU = 0.001;
            theta_p = ((U_main + dU) * seg_dist) / (2 * c_water);
            alpha_p = (flow_kg_s - theta_p) / (flow_kg_s + theta_p);
            T_in_perturbed = (T_prev - (1 - alpha_p) * T_soil) / alpha_p;
            dTdU = abs(T_in_perturbed - T_in) / dU;
            sigma_T_from_U = dTdU * U_main_uncertainty;

            sigma_T_in = sqrt(sigma_T_from_output^2 + sigma_T_from_U^2);

            T_junctions(sort_idx(c)) = T_in;
            T_uncertainties(sort_idx(c)) = sigma_T_in;

            T_prev = T_in;
            sigma_T_prev = sigma_T_in;
        end
    end
end