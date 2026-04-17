function save_daily_diagnostics(logger, ground_truth_csac, csac, output_folder, filter_name)
%SAVE_DAILY_DIAGNOSTICS Saves a compact daily-resolution CSV for analysis.
%   One row per house per day. Contains estimation errors, uncertainties, and NIS.
%   Only days with at least one filter update are included.

    if nargin < 5
        filter_name = 'ukf';
    end

    timestamps = logger.timestamps(:);
    num_houses = size(logger.state_estimates, 2);

    if isempty(timestamps) || all(isnat(timestamps))
        return;
    end

    % Extract series [time x houses]
    offset_est = squeeze(logger.state_estimates(1,:,:))';
    U_est      = squeeze(logger.state_estimates(2,:,:))';
    P_offset   = squeeze(logger.covariance_posterior(1,:,:))';  % variance (P diagonal)
    P_U        = squeeze(logger.covariance_posterior(2,:,:))';
    innovations = squeeze(logger.innovations)';                 % [time x houses]
    innov_var   = squeeze(logger.innovation_variance)';
    nis_values  = squeeze(logger.nis)';

    % Build day labels
    day_label = dateshift(timestamps, 'start', 'day');
    unique_days = unique(day_label(~isnat(day_label)));

    % Collect rows
    all_rows = {};

    for i = 1:num_houses
        house_id = logger.house_ids(i);
        gt_row = ground_truth_csac(ground_truth_csac.house_id == house_id, :);
        if isempty(gt_row), continue; end

        true_offset = gt_row.true_offset(1);
        true_U = gt_row.true_U(1);

        for d = 1:length(unique_days)
            mask = (day_label == unique_days(d)) & isfinite(offset_est(:,i));
            if ~any(mask), continue; end

            idx = find(mask);
            last_idx = idx(end);

            % Daily innovation stats (only finite values)
            day_innov = innovations(idx, i);
            day_innov = day_innov(isfinite(day_innov));

            day_nis = nis_values(idx, i);
            day_nis = day_nis(isfinite(day_nis));

            % Skip days with no actual filter updates
            if isempty(day_innov)
                continue;
            end

            row = { ...
                unique_days(d), ...
                house_id, ...
                round(offset_est(last_idx, i), 4), ...
                round(U_est(last_idx, i), 6), ...
                round(offset_est(last_idx, i) - true_offset, 4), ...
                round(U_est(last_idx, i) - true_U, 6), ...
                round(sqrt(max(0, P_offset(last_idx, i))), 4), ...
                round(sqrt(max(0, P_U(last_idx, i))), 6), ...
                round(mean(abs(day_innov)), 4), ...
                round(max(abs(day_innov)), 4), ...
                numel(day_innov), ...
                round(mean(day_nis), 4), ...
            };

            all_rows = [all_rows; row]; %#ok<AGROW>
        end
    end

    if isempty(all_rows)
        return;
    end

    T = cell2table(all_rows, 'VariableNames', { ...
        'date', 'house_id', ...
        'offset_est', 'U_est', ...
        'offset_err', 'U_err', ...
        'std_offset', 'std_U', ...
        'mean_abs_innov', 'max_abs_innov', 'n_updates', ...
        'mean_nis'});

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    filename = fullfile(output_folder, sprintf('%s_csac_%d_daily.csv', filter_name, csac));
    writetable(T, filename);
    fprintf('Daily diagnostics saved: %s\n', filename);
end