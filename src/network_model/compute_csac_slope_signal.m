function signal = compute_csac_slope_signal(cs, config)
%COMPUTE_CSAC_SLOPE_SIGNAL Computes position-dependent slope of offset and
%   U_service within a single CSAC.
%
%   Returns a signal struct with slope information that can be aggregated
%   across multiple CSACs for a shared U_csac estimate.
%
%   The offsets are de-meaned within the CSAC so that the overall offset
%   level does not affect the slope calculation.
%
%   Args:
%       cs (struct): CSAC state struct.
%       config (struct): Project configuration.
%
%   Returns:
%       signal (struct): Slope signal with fields:
%           .valid (logical): Whether enough data to compute slopes.
%           .num_valid_houses (int): Number of converged houses used.
%           .slope_offset (scalar): Slope of de-meaned offset vs position [°C/m].
%           .slope_U (scalar): Slope of U_service vs position [W/m/K/m].
%           .corr_offset (scalar): Correlation of de-meaned offset vs position.
%           .corr_U (scalar): Correlation of U_service vs position.
%           .pos_range (scalar): Range of house positions [m].
%           .total_gradient_offset (scalar): slope_offset * pos_range [°C].
%           .total_gradient_U (scalar): slope_U * pos_range [W/m/K].

    cfg = config.project.csac_U_estimation;
    num_houses = cs.num_houses;

    %% Extract estimates and positions
    offsets = zeros(num_houses, 1);
    U_services = zeros(num_houses, 1);
    positions = zeros(num_houses, 1);
    P_offsets = zeros(num_houses, 1);
    P_U = zeros(num_houses, 1);

    for i = 1:num_houses
        offsets(i) = cs.ukf_states{i}.x(1);
        U_services(i) = cs.ukf_states{i}.x(2);
        P_offsets(i) = cs.ukf_states{i}.P(1,1);
        P_U(i) = cs.ukf_states{i}.P(2,2);
        house_id = cs.house_ids(i);
        house_rows = cs.meter_data(cs.meter_data.house_id == house_id, :);
        if ~isempty(house_rows)
            positions(i) = house_rows.x_pos_m(1);
        end
    end

    %% Determine which houses have converged enough to use
    P_off_thresh = cfg.P_offset_convergence_threshold^2;
    P_U_thresh = cfg.P_U_convergence_threshold^2;

    valid_offset = P_offsets < P_off_thresh & P_offsets > 0;
    valid_U = P_U < P_U_thresh & P_U > 0;

    %% Check minimum houses
    if sum(valid_offset) < 3
        signal = make_empty_signal();
        return;
    end

    %% Compute offset slope (using de-meaned offsets)
    % De-mean removes the CSAC-level shared offset ambiguity
    w_off = zeros(num_houses, 1);
    w_off(valid_offset) = 1 ./ sqrt(P_offsets(valid_offset));
    w_off = w_off / sum(w_off);

    offsets_demeaned = offsets - sum(w_off .* offsets);
    [slope_offset, corr_offset] = weighted_linear_slope(positions, offsets_demeaned, w_off);

    %% Compute U_service slope
    if sum(valid_U) >= 3
        w_U = zeros(num_houses, 1);
        w_U(valid_U) = 1 ./ sqrt(P_U(valid_U));
        w_U = w_U / sum(w_U);
        [slope_U, corr_U] = weighted_linear_slope(positions, U_services, w_U);
    else
        slope_U = 0;
        corr_U = 0;
    end

    %% Compute total gradients
    pos_range = max(positions) - min(positions);

    %% Pack signal
    signal.valid = true;
    signal.num_valid_houses = sum(valid_offset);
    signal.slope_offset = slope_offset;
    signal.slope_U = slope_U;
    signal.corr_offset = corr_offset;
    signal.corr_U = corr_U;
    signal.pos_range = pos_range;
    signal.total_gradient_offset = slope_offset * pos_range;
    signal.total_gradient_U = slope_U * pos_range;
end


function [slope, correlation] = weighted_linear_slope(x, y, w)
%WEIGHTED_LINEAR_SLOPE Weighted linear regression slope and correlation.
    x_mean = sum(w .* x);
    y_mean = sum(w .* y);
    x_c = x - x_mean;
    y_c = y - y_mean;
    num = sum(w .* x_c .* y_c);
    den = sum(w .* x_c.^2);

    if abs(den) < 1e-12
        slope = 0;
        correlation = 0;
        return;
    end

    slope = num / den;
    std_x = sqrt(sum(w .* x_c.^2));
    std_y = sqrt(sum(w .* y_c.^2));
    if std_x > 0 && std_y > 0
        correlation = num / (std_x * std_y);
    else
        correlation = 0;
    end
end


function signal = make_empty_signal()
    signal.valid = false;
    signal.num_valid_houses = 0;
    signal.slope_offset = NaN;
    signal.slope_U = NaN;
    signal.corr_offset = NaN;
    signal.corr_U = NaN;
    signal.pos_range = NaN;
    signal.total_gradient_offset = NaN;
    signal.total_gradient_U = NaN;
end