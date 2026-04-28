function [U_csac_new, diagnostics] = estimate_csac_U(cs, config)
%ESTIMATE_CSAC_U Estimates the CSAC pipe U-value from position-offset correlation.
%
%   If the CSAC pipe U-value is wrong, the per-house offset estimates will
%   show a systematic trend with position along the CSAC:
%   - U_csac too high → offsets increase with distance from inlet
%   - U_csac too low  → offsets decrease with distance from inlet
%
%   This function detects that trend and adjusts U_csac accordingly.
%
%   Args:
%       cs (struct): CSAC state struct with current estimates.
%       config (struct): Project configuration.
%
%   Returns:
%       U_csac_new (scalar): Updated CSAC pipe U-value.
%       diagnostics (struct): Diagnostic information.

    %% Extract current per-house estimates and positions
    num_houses = cs.num_houses;
    offsets = zeros(num_houses, 1);
    positions = zeros(num_houses, 1);
    weights = zeros(num_houses, 1);
    P_offsets = zeros(num_houses, 1);

    for i = 1:num_houses
        offsets(i) = cs.ukf_states{i}.x(1);
        P_offsets(i) = cs.ukf_states{i}.P(1,1);
    end

    % Get positions from meter_data (x_pos_m)
    % We need unique positions per house, sorted by house_id
    for i = 1:num_houses
        house_id = cs.house_ids(i);
        house_rows = cs.meter_data(cs.meter_data.house_id == house_id, :);
        if ~isempty(house_rows)
            positions(i) = house_rows.x_pos_m(1);
        end
    end

    % Weight by inverse uncertainty — confident estimates matter more
    for i = 1:num_houses
        if P_offsets(i) > 0
            weights(i) = 1 / sqrt(P_offsets(i));
        else
            weights(i) = 0;
        end
    end

    % Normalize weights
    if sum(weights) > 0
        weights = weights / sum(weights);
    else
        U_csac_new = cs.U_csac;
        diagnostics.slope = NaN;
        diagnostics.correlation = NaN;
        diagnostics.adjustment = 0;
        return;
    end

    %% Compute weighted linear regression: offset = a + b * position
    % b > 0 means offsets increase with distance → U_csac too high
    % b < 0 means offsets decrease with distance → U_csac too low
    pos_mean = sum(weights .* positions);
    off_mean = sum(weights .* offsets);

    pos_centered = positions - pos_mean;
    off_centered = offsets - off_mean;

    numerator = sum(weights .* pos_centered .* off_centered);
    denominator = sum(weights .* pos_centered.^2);

    if abs(denominator) < 1e-12
        % All houses at same position or zero variance — can't estimate
        U_csac_new = cs.U_csac;
        diagnostics.slope = NaN;
        diagnostics.correlation = NaN;
        diagnostics.adjustment = 0;
        return;
    end

    slope = numerator / denominator;  % °C per meter

    % Correlation coefficient for diagnostics
    std_pos = sqrt(sum(weights .* pos_centered.^2));
    std_off = sqrt(sum(weights .* off_centered.^2));
    if std_pos > 0 && std_off > 0
        correlation = numerator / (std_pos * std_off);
    else
        correlation = 0;
    end

    %% Convert slope to U_csac adjustment
    % The relationship between slope and U_csac error depends on flow and
    % temperature conditions. We use a conservative proportional correction.
    %
    % Physical reasoning:
    % A slope of S [°C/m] in offset means the model is over/under-predicting
    % temperature drop by S [°C/m] along the CSAC pipe.
    % The temperature drop per meter is approximately:
    %   dT/dx ≈ U_csac / (flow_total * c) * (T - T_soil)
    % So a change in U_csac of dU causes a slope change of approximately:
    %   dS ≈ dU / (flow_total * c) * (T - T_soil)
    %
    % Rather than computing this exactly, we use a gain parameter.

    U_csac_gain = config.project.csac_U_estimation.gain;
    max_adjustment = config.project.csac_U_estimation.max_adjustment_per_step;
    U_min = config.project.csac_U_estimation.U_min;
    U_max = config.project.csac_U_estimation.U_max;
    min_correlation = config.project.csac_U_estimation.min_correlation;

    % Only adjust if the correlation is strong enough
    if abs(correlation) < min_correlation
        adjustment = 0;
    else
        % Positive slope → offsets increase with distance → U_csac too high
        % So we should DECREASE U_csac → adjustment is negative
        raw_adjustment = -U_csac_gain * slope;
        adjustment = max(-max_adjustment, min(max_adjustment, raw_adjustment));
    end

    U_csac_new = cs.U_csac + adjustment;
    U_csac_new = max(U_min, min(U_max, U_csac_new));

    %% Diagnostics
    diagnostics.slope = slope;
    diagnostics.correlation = correlation;
    diagnostics.adjustment = adjustment;
    diagnostics.offsets = offsets;
    diagnostics.positions = positions;
    diagnostics.weights = weights;
end