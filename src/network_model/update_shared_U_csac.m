function [U_csac_new, diagnostics] = update_shared_U_csac(all_cs, U_csac_current, config)
%UPDATE_SHARED_U_CSAC Updates the shared CSAC pipe U-value using slope
%   signals aggregated across all CSACs in the network.
%
%   Collects position-dependent slope signals from each CSAC, combines
%   them weighted by number of valid houses, and adjusts U_csac.
%
%   Args:
%       all_cs (cell array): CSAC state structs.
%       U_csac_current (scalar): Current shared U_csac value.
%       config (struct): Project configuration.
%
%   Returns:
%       U_csac_new (scalar): Updated U_csac value.
%       diagnostics (struct): Diagnostic information.

    cfg = config.project.csac_U_estimation;
    num_csacs = numel(all_cs);

    %% Collect slope signals from all CSACs
    signals = cell(num_csacs, 1);
    total_valid_houses = 0;
    weighted_gradient_offset = 0;
    weighted_gradient_U = 0;
    weighted_corr_offset = 0;

    for c = 1:num_csacs
        signals{c} = compute_csac_slope_signal(all_cs{c}, config);

        if signals{c}.valid
            n = signals{c}.num_valid_houses;
            total_valid_houses = total_valid_houses + n;
            weighted_gradient_offset = weighted_gradient_offset + n * signals{c}.total_gradient_offset;
            weighted_gradient_U = weighted_gradient_U + n * signals{c}.total_gradient_U;
            weighted_corr_offset = weighted_corr_offset + n * abs(signals{c}.corr_offset);
        end
    end

    %% Check if we have enough data
    if total_valid_houses < 6
        U_csac_new = U_csac_current;
        diagnostics.adjusted = false;
        diagnostics.reason = 'insufficient_houses';
        diagnostics.total_valid_houses = total_valid_houses;
        diagnostics.signals = signals;
        return;
    end

    %% Compute aggregate signals
    avg_gradient_offset = weighted_gradient_offset / total_valid_houses;
    avg_gradient_U = weighted_gradient_U / total_valid_houses;
    avg_corr = weighted_corr_offset / total_valid_houses;

    %% Check if signal is strong enough
    if avg_corr < cfg.min_correlation
        U_csac_new = U_csac_current;
        diagnostics.adjusted = false;
        diagnostics.reason = 'weak_correlation';
        diagnostics.avg_corr = avg_corr;
        diagnostics.avg_gradient_offset = avg_gradient_offset;
        diagnostics.total_valid_houses = total_valid_houses;
        diagnostics.signals = signals;
        return;
    end

    %% Compute adjustment
    % Positive offset gradient → U_csac too high → decrease
    raw_adjustment = cfg.gain * avg_gradient_offset;

    % Boost if U_service slope confirms (opposite sign to offset slope)
    if sign(avg_gradient_U) == sign(avg_gradient_offset) && abs(avg_gradient_U) > 0.0001
        raw_adjustment = raw_adjustment * 1.5;
    end

    % Clamp
    adjustment = max(-cfg.max_adjustment_per_step, min(cfg.max_adjustment_per_step, raw_adjustment));

    % Apply
    U_csac_new = U_csac_current + adjustment;
    U_csac_new = max(cfg.U_min, min(cfg.U_max, U_csac_new));

    %% Diagnostics
    diagnostics.adjusted = true;
    diagnostics.reason = 'updated';
    diagnostics.avg_gradient_offset = avg_gradient_offset;
    diagnostics.avg_gradient_U = avg_gradient_U;
    diagnostics.avg_corr = avg_corr;
    diagnostics.raw_adjustment = raw_adjustment;
    diagnostics.adjustment = adjustment;
    diagnostics.total_valid_houses = total_valid_houses;
    diagnostics.U_old = U_csac_current;
    diagnostics.U_new = U_csac_new;
    diagnostics.signals = signals;
end