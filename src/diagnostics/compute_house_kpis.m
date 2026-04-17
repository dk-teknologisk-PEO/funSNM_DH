function kpis = compute_house_kpis(state_estimates, covariance_posterior, timestamps, ...
    true_offset_traj, true_U_traj, kpi_config)
%COMPUTE_HOUSE_KPIS Computes performance KPIs for a single house.
%
%   Args:
%       state_estimates (2xT matrix): Estimated [offset; U] at each timestep.
%       covariance_posterior (2xT matrix): Posterior variances [P_offset; P_U] 
%           at each timestep (diagonal elements of P).
%       timestamps (1xT datetime): Timestamp for each step.
%       true_offset_traj (1xT double): True offset at each timestep.
%       true_U_traj (1xT double): True U-value at each timestep.
%       kpi_config (struct): KPI parameters:
%           .offset_tolerance (scalar): Error band for offset [°C].
%           .U_tolerance (scalar): Error band for U [W/m/K].
%           .convergence_P_offset (scalar): P threshold for offset convergence [°C].
%           .convergence_P_U (scalar): P threshold for U convergence [W/m/K].
%           .convergence_hold_days (int): Active days P must stay below threshold.
%
%   Returns:
%       kpis (struct): Computed KPIs.

    T = length(timestamps);

    %% Compute error and uncertainty trajectories
    err_offset = state_estimates(1, :) - true_offset_traj;
    err_U = state_estimates(2, :) - true_U_traj;
    std_offset = sqrt(covariance_posterior(1, :));
    std_U = sqrt(covariance_posterior(2, :));

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

    %% KPI 2: Convergence Time (uncertainty-based)
    % Convergence = when the posterior std drops below threshold and stays there
    hold_hours = kpi_config.convergence_hold_days * 24;

    [kpis.convergence_days_offset, conv_idx_offset] = find_convergence_time_P(...
        std_offset, dt_hours, timestamps, kpi_config.convergence_P_offset, ...
        hold_hours, active_idx);

    [kpis.convergence_days_U, conv_idx_U] = find_convergence_time_P(...
        std_U, dt_hours, timestamps, kpi_config.convergence_P_U, ...
        hold_hours, active_idx);

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

    %% KPI 6: End-of-Season Errors and Uncertainties
    kpis.end_of_season_errors = compute_end_of_season_errors(...
        state_estimates, covariance_posterior, err_offset, err_U, ...
        is_active, timestamps);

    %% KPI 7: Initial and final state snapshots (for bar plots)
    first_idx = active_idx(1);
    last_idx = active_idx(end);
    kpis.initial_state.offset = state_estimates(1, first_idx);
    kpis.initial_state.U = state_estimates(2, first_idx);
    kpis.initial_state.std_offset = std_offset(first_idx);
    kpis.initial_state.std_U = std_U(first_idx);
    kpis.final_state.offset = state_estimates(1, last_idx);
    kpis.final_state.U = state_estimates(2, last_idx);
    kpis.final_state.std_offset = std_offset(last_idx);
    kpis.final_state.std_U = std_U(last_idx);

    %% KPI 8: Final absolute error
    kpis.final_error_offset = err_offset(last_idx);
    kpis.final_error_U = err_U(last_idx);
end


%% ===== LOCAL HELPER FUNCTIONS =====

function [conv_days, conv_idx] = find_convergence_time_P(std_traj, dt_hours, timestamps, ...
    threshold, hold_hours, active_idx)
%FIND_CONVERGENCE_TIME_P Finds when posterior std first drops below threshold and stays.
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

        if std_traj(si) > threshold
            continue;
        end

        % Walk forward accumulating active time below threshold
        accumulated_hours = 0;
        broke_out = false;

        for check_k = (start_k + 1):n_active
            ci = active_idx(check_k);

            if std_traj(ci) > threshold
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

        % Reached end of data while still below threshold
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


function eos = compute_end_of_season_errors(state_est, cov_post, err_offset, err_U, ...
    is_active, timestamps)
%COMPUTE_END_OF_SEASON_ERRORS Records state, uncertainty and error at the
%   last active timestep before each gap of >30 days, plus final timestep.

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
    season_ends(end+1) = active_idx(end);

    n_seasons = numel(season_ends);
    season_num = (1:n_seasons)';
    end_dates = timestamps(season_ends)';
    offset_est = state_est(1, season_ends)';
    U_est = state_est(2, season_ends)';
    offset_std = sqrt(cov_post(1, season_ends))';
    U_std = sqrt(cov_post(2, season_ends))';
    offset_err = err_offset(season_ends)';
    U_err = err_U(season_ends)';

    eos = table(season_num, end_dates, offset_est, offset_std, offset_err, ...
        U_est, U_std, U_err, ...
        'VariableNames', {'season', 'end_date', 'offset_est', 'offset_std', ...
        'offset_error_C', 'U_est', 'U_std', 'U_error_WmK'});
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
    kpis.initial_state = struct('offset', NaN, 'U', NaN, 'std_offset', NaN, 'std_U', NaN);
    kpis.final_state = struct('offset', NaN, 'U', NaN, 'std_offset', NaN, 'std_U', NaN);
    kpis.final_error_offset = NaN;
    kpis.final_error_U = NaN;
end