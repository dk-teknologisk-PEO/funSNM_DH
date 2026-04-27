function kpis = compute_house_kpis(state_estimates, covariance_posterior, timestamps, ...
    true_offset_traj, true_U_traj, kpi_config)
%COMPUTE_HOUSE_KPIS Computes performance KPIs for a single house.

    % --- Force all 1D inputs to row vectors ---
    timestamps = timestamps(:)';
    true_offset_traj = true_offset_traj(:)';
    true_U_traj = true_U_traj(:)';

    if size(state_estimates,1) ~= 2 && size(state_estimates,2) == 2
        state_estimates = state_estimates';
    end
    if size(covariance_posterior,1) ~= 2 && size(covariance_posterior,2) == 2
        covariance_posterior = covariance_posterior';
    end

    T = length(timestamps);

    %% Compute error and uncertainty trajectories
    err_offset = state_estimates(1, :) - true_offset_traj;
    err_U = state_estimates(2, :) - true_U_traj;
    std_offset = sqrt(covariance_posterior(1, :));
    std_U = sqrt(covariance_posterior(2, :));

    %% Define valid timesteps
    % Valid = timestamp is not NaT AND state is not NaN
    is_valid_time = ~isnat(timestamps);
    is_valid_state = ~isnan(err_offset) & ~isnan(err_U);
    is_valid = is_valid_time & is_valid_state;
    valid_idx = find(is_valid);

    if numel(valid_idx) < 2
        kpis = empty_kpis();
        return;
    end

    %% Compute time deltas between valid timesteps (in hours)
    dt_hours = zeros(1, T);
    for k = 2:numel(valid_idx)
        idx = valid_idx(k);
        idx_prev = valid_idx(k-1);
        dt_hours(idx) = hours(timestamps(idx) - timestamps(idx_prev));
    end
    dt_hours(valid_idx(1)) = 1;

    % Cap large gaps so off-season periods do not dominate
    max_dt_hours = 48;
    dt_hours = min(dt_hours, max_dt_hours);

    total_time = sum(dt_hours(is_valid));

    %% KPI 1: Time-Weighted MAE
    if total_time > 0
        kpis.tw_mae_offset = sum(dt_hours(is_valid) .* abs(err_offset(is_valid))) / total_time;
        kpis.tw_mae_U = sum(dt_hours(is_valid) .* abs(err_U(is_valid))) / total_time;
    else
        kpis.tw_mae_offset = NaN;
        kpis.tw_mae_U = NaN;
    end

    %% KPI 2: Convergence Time (uncertainty-based)
    hold_hours = kpi_config.convergence_hold_days * 24;

    [kpis.convergence_days_offset, conv_idx_offset] = find_convergence_time_P(...
        std_offset, dt_hours, timestamps, kpi_config.convergence_P_offset, ...
        hold_hours, valid_idx);

    [kpis.convergence_days_U, conv_idx_U] = find_convergence_time_P(...
        std_U, dt_hours, timestamps, kpi_config.convergence_P_U, ...
        hold_hours, valid_idx);

    %% KPI 3: Steady-State MAE (after convergence)
    kpis.steady_state_mae_offset = compute_post_conv_weighted_mae(...
        err_offset, is_valid, dt_hours, conv_idx_offset, T);
    kpis.steady_state_mae_U = compute_post_conv_weighted_mae(...
        err_U, is_valid, dt_hours, conv_idx_U, T);

    %% KPI 4: Stability Index (std of error after convergence)
    kpis.stability_index_offset = compute_post_conv_std(err_offset, is_valid, conv_idx_offset, T);
    kpis.stability_index_U = compute_post_conv_std(err_U, is_valid, conv_idx_U, T);

    %% KPI 5: Maximum Excursion After Convergence
    kpis.max_excursion_offset = compute_post_conv_max(err_offset, is_valid, conv_idx_offset, T);
    kpis.max_excursion_U = compute_post_conv_max(err_U, is_valid, conv_idx_U, T);

    %% KPI 6: End-of-Season Errors and Uncertainties
    kpis.end_of_season_errors = compute_end_of_season_errors(...
        state_estimates, covariance_posterior, err_offset, err_U, ...
        is_valid, timestamps);

    %% KPI 7: Initial and final state snapshots
    first_idx = valid_idx(1);
    last_idx = valid_idx(end);

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
    threshold, hold_hours, valid_idx)

    conv_days = NaN;
    conv_idx = NaN;

    if numel(valid_idx) < 2
        return;
    end

    first_valid_time = timestamps(valid_idx(1));

    for start_k = 1:numel(valid_idx)
        si = valid_idx(start_k);

        if std_traj(si) > threshold
            continue;
        end

        accumulated_hours = 0;
        broke_out = false;

        for check_k = (start_k + 1):numel(valid_idx)
            ci = valid_idx(check_k);

            if std_traj(ci) > threshold
                broke_out = true;
                break;
            end

            accumulated_hours = accumulated_hours + dt_hours(ci);

            if accumulated_hours >= hold_hours
                conv_idx = si;
                conv_days = days(timestamps(si) - first_valid_time);
                return;
            end
        end

        if ~broke_out && accumulated_hours >= hold_hours
            conv_idx = si;
            conv_days = days(timestamps(si) - first_valid_time);
            return;
        end
    end
end


function val = compute_post_conv_weighted_mae(err, is_valid, dt_hours, conv_idx, T)
    if ~isfinite(conv_idx)
        val = NaN;
        return;
    end

    err = err(:)';
    is_valid = is_valid(:)';
    dt_hours = dt_hours(:)';

    idx_vec = 1:T;
    post_conv = is_valid & (idx_vec >= conv_idx);

    total_time = sum(dt_hours(post_conv));

    if total_time <= 0
        val = NaN;
    else
        val = sum(dt_hours(post_conv) .* abs(err(post_conv))) / total_time;
    end
end


function val = compute_post_conv_std(err, is_valid, conv_idx, T)
    if ~isfinite(conv_idx)
        val = NaN;
        return;
    end

    err = err(:)';
    is_valid = is_valid(:)';

    idx_vec = 1:T;
    post_conv = is_valid & (idx_vec >= conv_idx);

    if sum(post_conv) < 2
        val = NaN;
    else
        val = std(err(post_conv));
    end
end


function val = compute_post_conv_max(err, is_valid, conv_idx, T)
    if ~isfinite(conv_idx)
        val = NaN;
        return;
    end

    err = err(:)';
    is_valid = is_valid(:)';

    idx_vec = 1:T;
    post_conv = is_valid & (idx_vec >= conv_idx);

    if sum(post_conv) == 0
        val = NaN;
    else
        val = max(abs(err(post_conv)));
    end
end


function eos = compute_end_of_season_errors(state_est, cov_post, err_offset, err_U, ...
    is_valid, timestamps)

    timestamps = timestamps(:)';
    is_valid = is_valid(:)';
    err_offset = err_offset(:)';
    err_U = err_U(:)';

    valid_idx = find(is_valid);
    if numel(valid_idx) < 2
        eos = table();
        return;
    end

    gap_threshold_days = 30;
    season_ends = [];

    for k = 1:(numel(valid_idx) - 1)
        gap = days(timestamps(valid_idx(k+1)) - timestamps(valid_idx(k)));
        if gap > gap_threshold_days
            season_ends(end+1) = valid_idx(k); %#ok<AGROW>
        end
    end
    season_ends(end+1) = valid_idx(end);

    season_num = (1:numel(season_ends))';
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