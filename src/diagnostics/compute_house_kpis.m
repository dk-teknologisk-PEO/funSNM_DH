function kpis = compute_house_kpis(state_estimates, timestamps, true_offset_traj, true_U_traj, kpi_config)
%COMPUTE_HOUSE_KPIS Computes performance KPIs for a single house.
%
%   Args:
%       state_estimates (2xT matrix): Estimated [offset; U] at each timestep.
%       timestamps (1xT datetime): Timestamp for each step.
%       true_offset_traj (1xT double): True offset at each timestep.
%       true_U_traj (1xT double): True U-value at each timestep.
%       kpi_config (struct): KPI parameters:
%           .offset_tolerance (scalar): Convergence band for offset [°C].
%           .U_tolerance (scalar): Convergence band for U [W/m/K].
%           .convergence_hold_days (int): Active days estimate must stay in band.
%
%   Returns:
%       kpis (struct): Computed KPIs.

    T = length(timestamps);

    %% Compute error trajectories
    err_offset = state_estimates(1, :) - true_offset_traj;
    err_U = state_estimates(2, :) - true_U_traj;

    %% Find active timesteps
    % Active = timestamp is valid AND state has been updated at least once
    is_valid = ~isnat(timestamps);
    state_changed = [false, any(diff(state_estimates, 1, 2) ~= 0, 1)];
    has_been_updated = cumsum(state_changed) > 0;
    is_active = is_valid & has_been_updated;

    active_idx = find(is_active);
    if numel(active_idx) < 2
        kpis = empty_kpis();
        return;
    end

    %% Compute time deltas between active timesteps (in hours)
    dt_hours = zeros(1, T);
    for k = 2:numel(active_idx)
        idx = active_idx(k);
        idx_prev = active_idx(k-1);
        dt_hours(idx) = hours(timestamps(idx) - timestamps(idx_prev));
    end
    dt_hours(active_idx(1)) = 1; % nominal weight for first active step

    total_time = sum(dt_hours(is_active));

    %% KPI 1: Time-Weighted MAE
    kpis.tw_mae_offset = sum(dt_hours(is_active) .* abs(err_offset(is_active))) / total_time;
    kpis.tw_mae_U = sum(dt_hours(is_active) .* abs(err_U(is_active))) / total_time;

    %% KPI 2: Convergence Time
    hold_hours = kpi_config.convergence_hold_days * 24;

    [kpis.convergence_days_offset, conv_idx_offset] = find_convergence_time(...
        err_offset, dt_hours, timestamps, kpi_config.offset_tolerance, hold_hours, active_idx);

    [kpis.convergence_days_U, conv_idx_U] = find_convergence_time(...
        err_U, dt_hours, timestamps, kpi_config.U_tolerance, hold_hours, active_idx);

    %% KPI 3: Steady-State MAE (after convergence)
    kpis.steady_state_mae_offset = compute_post_conv_weighted_mae(...
        err_offset, is_active, dt_hours, conv_idx_offset, T);
    kpis.steady_state_mae_U = compute_post_conv_weighted_mae(...
        err_U, is_active, dt_hours, conv_idx_U, T);

    %% KPI 4: Stability Index (std of error after convergence)
    kpis.stability_index_offset = compute_post_conv_std(err_offset, is_active, conv_idx_offset, T);
    kpis.stability_index_U = compute_post_conv_std(err_U, is_active, conv_idx_U, T);

    %% KPI 5: Maximum Excursion After Convergence
    kpis.max_excursion_offset = compute_post_conv_max(err_offset, is_active, conv_idx_offset, T);
    kpis.max_excursion_U = compute_post_conv_max(err_U, is_active, conv_idx_U, T);

    %% KPI 6: End-of-Season Errors
    kpis.end_of_season_errors = compute_end_of_season_errors(...
        err_offset, err_U, is_active, timestamps);
end


%% ===== LOCAL HELPER FUNCTIONS =====

function [conv_days, conv_idx] = find_convergence_time(err, dt_hours, timestamps, tolerance, hold_hours, active_idx)
%FIND_CONVERGENCE_TIME Finds when the error first stays within tolerance.
%   Returns days from first active timestep to convergence, and the
%   timestep index where convergence occurs. NaN if never converges.

    conv_days = NaN;
    conv_idx = NaN;
    n_active = numel(active_idx);

    if n_active < 2
        return;
    end

    first_active_time = timestamps(active_idx(1));

    for start_k = 1:n_active
        si = active_idx(start_k);

        if abs(err(si)) > tolerance
            continue;
        end

        % Walk forward accumulating active time within tolerance
        accumulated_hours = 0;
        broke_out = false;

        for check_k = (start_k + 1):n_active
            ci = active_idx(check_k);

            if abs(err(ci)) > tolerance
                broke_out = true;
                break;
            end

            accumulated_hours = accumulated_hours + dt_hours(ci);

            if accumulated_hours >= hold_hours
                conv_idx = si;
                conv_days = days(timestamps(si) - first_active_time);
                return;
            end
        end

        % Reached end of data while still in tolerance
        if ~broke_out && accumulated_hours >= hold_hours
            conv_idx = si;
            conv_days = days(timestamps(si) - first_active_time);
            return;
        end
    end
end


function val = compute_post_conv_weighted_mae(err, is_active, dt_hours, conv_idx, T)
%COMPUTE_POST_CONV_WEIGHTED_MAE Time-weighted MAE after convergence.
    if ~isfinite(conv_idx)
        val = NaN;
        return;
    end
    post_conv = is_active & (1:T) >= conv_idx;
    total_time = sum(dt_hours(post_conv));
    if total_time == 0
        val = NaN;
    else
        val = sum(dt_hours(post_conv) .* abs(err(post_conv))) / total_time;
    end
end


function val = compute_post_conv_std(err, is_active, conv_idx, T)
%COMPUTE_POST_CONV_STD Standard deviation of error after convergence.
    if ~isfinite(conv_idx)
        val = NaN;
        return;
    end
    post_conv = is_active & (1:T) >= conv_idx;
    if sum(post_conv) < 2
        val = NaN;
    else
        val = std(err(post_conv));
    end
end


function val = compute_post_conv_max(err, is_active, conv_idx, T)
%COMPUTE_POST_CONV_MAX Maximum absolute error after convergence.
    if ~isfinite(conv_idx)
        val = NaN;
        return;
    end
    post_conv = is_active & (1:T) >= conv_idx;
    if sum(post_conv) == 0
        val = NaN;
    else
        val = max(abs(err(post_conv)));
    end
end


function eos = compute_end_of_season_errors(err_offset, err_U, is_active, timestamps)
%COMPUTE_END_OF_SEASON_ERRORS Records error at the last active timestep
%   before each gap of >30 days, plus the final active timestep.

    active_idx = find(is_active);
    if numel(active_idx) < 2
        eos = table();
        return;
    end

    gap_threshold_days = 30;
    season_ends = [];

    for k = 1:(numel(active_idx) - 1)
        gap = days(timestamps(active_idx(k+1)) - timestamps(active_idx(k)));
        if gap > gap_threshold_days
            season_ends(end+1) = active_idx(k); %#ok<AGROW>
        end
    end
    % Last active timestep is always a season end
    season_ends(end+1) = active_idx(end);

    n_seasons = numel(season_ends);
    season_num = (1:n_seasons)';
    end_dates = timestamps(season_ends)';
    offset_err = err_offset(season_ends)';
    U_err = err_U(season_ends)';

    eos = table(season_num, end_dates, offset_err, U_err, ...
        'VariableNames', {'season', 'end_date', 'offset_error_C', 'U_error_WmK'});
end


function kpis = empty_kpis()
%EMPTY_KPIS Returns a KPI struct with all NaN values for insufficient data.
    kpis.tw_mae_offset = NaN;
    kpis.tw_mae_U = NaN;
    kpis.convergence_days_offset = NaN;
    kpis.convergence_days_U = NaN;
    kpis.steady_state_mae_offset = NaN;
    kpis.steady_state_mae_U = NaN;
    kpis.stability_index_offset = NaN;
    kpis.stability_index_U = NaN;
    kpis.max_excursion_offset = NaN;
    kpis.max_excursion_U = NaN;
    kpis.end_of_season_errors = table();
end