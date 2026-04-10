function save_diagnostic_summary_detailed(logger, ground_truth_csac, csac, output_folder, filter_name)
%SAVE_DIAGNOSTIC_SUMMARY_DETAILED
% Creates a per-house diagnostic CSV with focus on:
%   1) how good winter performance is,
%   2) when the estimates first become good,
%   3) how long the good/stable period lasts,
%   4) when divergence starts.
%
% MATLAB R2023 compatible.

    if nargin < 5
        filter_name = 'ukf';
    end

    timestamps = logger.timestamps(:);
    num_houses = size(logger.state_estimates, 2);
    num_steps = size(logger.state_estimates, 3);

    if isempty(timestamps) || all(isnat(timestamps))
        warning('No valid timestamps found for CSAC %d. Detailed summary not saved.', csac);
        return;
    end

    % --- Adjustable thresholds for "good performance" ---
    offset_good_thresh = 0.25;   % degC
    U_good_thresh      = 0.015;  % W/m/K

    % Divergence thresholds
    offset_div_thresh = 0.75;    % degC
    U_div_thresh      = 0.03;    % W/m/K

    % Require sustained good/bad behavior over a small window
    sustain_len = 5;

    % Winter period definition
    is_winter = month(timestamps) <= 3;   % Jan-Mar inclusive
    is_feb_to_mar = (month(timestamps) >= 2) & (month(timestamps) <= 3);

    % Tail period: last 15% of valid timesteps
    tail_frac = 0.15;

    % State bounds inferred from actual estimates if config not available here
    offset_series_all = squeeze(logger.state_estimates(1,:,:)); % houses x time after squeeze may vary
    U_series_all      = squeeze(logger.state_estimates(2,:,:));

    % Ensure dimensions are [num_houses x num_steps]
    if size(offset_series_all,1) ~= num_houses
        offset_series_all = offset_series_all';
    end
    if size(U_series_all,1) ~= num_houses
        U_series_all = U_series_all';
    end

    summary = table();

    for i = 1:num_houses
        house_id = logger.house_ids(i);

        gt_row = ground_truth_csac(ground_truth_csac.house_id == house_id, :);
        if isempty(gt_row)
            continue
        end

        true_offset = gt_row.true_offset(1);
        true_U = gt_row.true_U(1);

        offset_est = offset_series_all(i, :)';
        U_est      = U_series_all(i, :)';

        valid = isfinite(offset_est) & isfinite(U_est) & ~isnat(timestamps);

        if ~any(valid)
            continue
        end

        t_valid = timestamps(valid);
        offset_est_valid = offset_est(valid);
        U_est_valid      = U_est(valid);

        offset_err = offset_est_valid - true_offset;
        U_err      = U_est_valid - true_U;

        % --- Final values ---
        final_offset = offset_est_valid(end);
        final_U      = U_est_valid(end);
        final_offset_error = offset_err(end);
        final_U_error      = U_err(end);

        % --- Overall RMSE ---
        offset_rmse_all = sqrt(mean(offset_err.^2, 'omitnan'));
        U_rmse_all      = sqrt(mean(U_err.^2, 'omitnan'));

        % --- Winter RMSE ---
        valid_winter = valid & is_winter;
        if any(valid_winter)
            ow = offset_est(valid_winter) - true_offset;
            uw = U_est(valid_winter) - true_U;
            offset_rmse_winter = sqrt(mean(ow.^2, 'omitnan'));
            U_rmse_winter      = sqrt(mean(uw.^2, 'omitnan'));
        else
            offset_rmse_winter = NaN;
            U_rmse_winter      = NaN;
        end

        % --- Feb-Mar RMSE ---
        valid_feb_to_mar = valid & is_feb_to_mar;
        if any(valid_feb_to_mar)
            ofm = offset_est(valid_feb_to_mar) - true_offset;
            ufm = U_est(valid_feb_to_mar) - true_U;
            offset_rmse_feb_to_mar = sqrt(mean(ofm.^2, 'omitnan'));
            U_rmse_feb_to_mar      = sqrt(mean(ufm.^2, 'omitnan'));
        else
            offset_rmse_feb_to_mar = NaN;
            U_rmse_feb_to_mar      = NaN;
        end

        % --- Tail metrics ---
        n_valid = numel(offset_est_valid);
        tail_start = max(1, floor((1-tail_frac)*n_valid));
        offset_tail = offset_est_valid(tail_start:end);
        U_tail      = U_est_valid(tail_start:end);

        offset_tail_std = std(offset_tail, 'omitnan');
        U_tail_std      = std(U_tail, 'omitnan');

        % --- Jump metrics ---
        doffset = abs(diff(offset_est_valid));
        dU      = abs(diff(U_est_valid));

        if isempty(doffset)
            max_offset_jump = NaN;
            jump_idx_offset = NaN;
        else
            [max_offset_jump, idx1] = max(doffset);
            jump_idx_offset = idx1 + 1;
        end

        if isempty(dU)
            max_U_jump = NaN;
            jump_idx_U = NaN;
        else
            [max_U_jump, idx2] = max(dU);
            jump_idx_U = idx2 + 1;
        end

        % --- Good period detection ---
        good_mask = (abs(offset_err) <= offset_good_thresh) & (abs(U_err) <= U_good_thresh);
        good_mask = sustained_mask(good_mask, sustain_len);

        [good_start_idx, good_end_idx, good_len] = longest_true_run(good_mask);

        if isnan(good_start_idx)
            good_start_time = NaT;
            good_end_time   = NaT;
        else
            good_start_time = t_valid(good_start_idx);
            good_end_time   = t_valid(good_end_idx);
        end

        % --- First divergence after first good period ---
        div_mask = (abs(offset_err) >= offset_div_thresh) | (abs(U_err) >= U_div_thresh);
        div_mask = sustained_mask(div_mask, sustain_len);

        divergence_idx = NaN;
        divergence_time = NaT;

        if ~isnan(good_end_idx) && good_end_idx < numel(div_mask)
            idx_rel = find(div_mask(good_end_idx+1:end), 1, 'first');
            if ~isempty(idx_rel)
                divergence_idx = good_end_idx + idx_rel;
                divergence_time = t_valid(divergence_idx);
            end
        end

        % --- Time spent "good" in winter ---
        winter_valid_local = month(t_valid) <= 3;
        if any(winter_valid_local)
            good_fraction_winter = mean(good_mask(winter_valid_local));
        else
            good_fraction_winter = NaN;
        end

        % --- Bias metrics in winter ---
        if any(winter_valid_local)
            mean_offset_error_winter = mean(offset_err(winter_valid_local), 'omitnan');
            mean_U_error_winter      = mean(U_err(winter_valid_local), 'omitnan');
        else
            mean_offset_error_winter = NaN;
            mean_U_error_winter      = NaN;
        end

        % --- Add row ---
        row = table( ...
            house_id, ...
            final_offset, final_U, true_offset, true_U, ...
            final_offset_error, final_U_error, ...
            offset_rmse_all, U_rmse_all, ...
            offset_rmse_winter, U_rmse_winter, ...
            offset_rmse_feb_to_mar, U_rmse_feb_to_mar, ...
            mean_offset_error_winter, mean_U_error_winter, ...
            max_offset_jump, jump_idx_offset, ...
            max_U_jump, jump_idx_U, ...
            offset_tail_std, U_tail_std, ...
            good_len, good_start_time, good_end_time, ...
            good_fraction_winter, ...
            divergence_idx, divergence_time, ...
            'VariableNames', { ...
            'house_id', ...
            'final_offset', 'final_U', 'true_offset', 'true_U', ...
            'offset_error_final', 'U_error_final', ...
            'offset_rmse_all', 'U_rmse_all', ...
            'offset_rmse_winter', 'U_rmse_winter', ...
            'offset_rmse_feb_to_mar', 'U_rmse_feb_to_mar', ...
            'mean_offset_error_winter', 'mean_U_error_winter', ...
            'max_offset_jump', 'jump_idx_offset', ...
            'max_U_jump', 'jump_idx_U', ...
            'offset_tail_std', 'U_tail_std', ...
            'good_period_len', 'good_period_start', 'good_period_end', ...
            'good_fraction_winter', ...
            'divergence_idx', 'divergence_time'});

        summary = [summary; row]; %#ok<AGROW>
    end

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    filename = fullfile(output_folder, sprintf('%s_csac_%d_summary_detailed.csv', filter_name, csac));
    writetable(summary, filename);

    fprintf('Detailed diagnostic summary saved: %s\n', filename);
end


function mask_out = sustained_mask(mask_in, sustain_len)
% Keep only regions where mask is true for at least sustain_len consecutive samples.
    mask_in = logical(mask_in(:));
    mask_out = false(size(mask_in));

    if isempty(mask_in)
        return
    end

    d = diff([false; mask_in; false]);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;

    run_lengths = run_ends - run_starts + 1;

    for k = 1:numel(run_starts)
        if run_lengths(k) >= sustain_len
            mask_out(run_starts(k):run_ends(k)) = true;
        end
    end
end


function [start_idx, end_idx, run_len] = longest_true_run(mask)
% Find longest consecutive run of true values.
    mask = logical(mask(:));

    start_idx = NaN;
    end_idx = NaN;
    run_len = 0;

    if isempty(mask) || ~any(mask)
        return
    end

    d = diff([false; mask; false]);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;
    run_lengths = run_ends - run_starts + 1;

    [run_len, idx] = max(run_lengths);
    start_idx = run_starts(idx);
    end_idx   = run_ends(idx);
end