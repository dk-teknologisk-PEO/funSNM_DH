function summary = compute_and_save_network_kpis(cs, csac_id, output_folder, kpi_config, true_trajectories)
%COMPUTE_AND_SAVE_NETWORK_KPIS Computes KPIs for all houses in a CSAC,
%   saves CSV results and generates visualization plots.
%
%   Args:
%       cs (struct): CSAC state struct after processing.
%       csac_id (int): CSAC identifier.
%       output_folder (char): Output directory for CSV and plot files.
%       kpi_config (struct): KPI parameters.
%       true_trajectories (cell, optional): Cell array of structs per house
%           with fields .offset (1xT) and .U (1xT). If not provided or
%           empty, constant trajectories are built from ground_truth.
%
%   Returns:
%       summary (table): Summary table with one row per house.

    num_houses = cs.num_houses;
    T = length(cs.logger.timestamps);

    % Build constant trajectories if not provided
    if nargin < 5 || isempty(true_trajectories)
        true_trajectories = cell(num_houses, 1);
        for i = 1:num_houses
            true_trajectories{i}.offset = repmat(cs.ground_truth.true_offset(i), 1, T);
            true_trajectories{i}.U = repmat(cs.ground_truth.true_U(i), 1, T);
        end
    end

    % Pre-allocate summary columns
    house_ids = cs.house_ids(:);
    tw_mae_offset = nan(num_houses, 1);
    tw_mae_U = nan(num_houses, 1);
    conv_days_offset = nan(num_houses, 1);
    conv_days_U = nan(num_houses, 1);
    ss_mae_offset = nan(num_houses, 1);
    ss_mae_U = nan(num_houses, 1);
    stability_offset = nan(num_houses, 1);
    stability_U = nan(num_houses, 1);
    max_exc_offset = nan(num_houses, 1);
    max_exc_U = nan(num_houses, 1);
    true_offset_initial = nan(num_houses, 1);
    true_offset_final = nan(num_houses, 1);
    true_offset_drift = nan(num_houses, 1);
    true_U_initial = nan(num_houses, 1);
    true_U_final = nan(num_houses, 1);
    final_err_offset = nan(num_houses, 1);
    final_err_U = nan(num_houses, 1);

    % For bar plot
    initial_offsets = nan(num_houses, 1);
    initial_U = nan(num_houses, 1);
    initial_std_offset = nan(num_houses, 1);
    initial_std_U = nan(num_houses, 1);
    final_offsets = nan(num_houses, 1);
    final_U_vals = nan(num_houses, 1);
    final_std_offset = nan(num_houses, 1);
    final_std_U = nan(num_houses, 1);

    all_eos = table();

    for i = 1:num_houses
        state_est = squeeze(cs.logger.state_estimates(:, i, :)); % 2xT
        cov_post = squeeze(cs.logger.covariance_posterior(:, i, :)); % 2xT

        kpis = compute_house_kpis(state_est, cov_post, cs.logger.timestamps, ...
            true_trajectories{i}.offset, true_trajectories{i}.U, kpi_config);

        % % TEMPORARY DEBUG — remove after fixing
        % if i == 1
        %     fprintf('  DEBUG house %d: T=%d, valid=%d, tw_mae=%.4f, conv_idx_off=%s, ss_mae=%.4f\n', ...
        %         house_ids(i), T, sum(~isnat(cs.logger.timestamps)), ...
        %         kpis.tw_mae_offset, ...
        %         mat2str(kpis.convergence_days_offset), ...
        %         kpis.steady_state_mae_offset);
        %     fprintf('  DEBUG: timestamps class=%s, size=%s\n', ...
        %         class(cs.logger.timestamps), mat2str(size(cs.logger.timestamps)));
        %     fprintf('  DEBUG: state_est size=%s, cov_post size=%s\n', ...
        %         mat2str(size(state_est)), mat2str(size(cov_post)));
        %     fprintf('  DEBUG: true_offset size=%s\n', ...
        %         mat2str(size(true_trajectories{i}.offset)));
        % end
        % if i == 1
        %     % Check for NaN in error trajectories
        %     err_check = state_est(1,:) - true_trajectories{i}.offset;
        %     fprintf('  DEBUG: err_offset has %d NaN out of %d\n', sum(isnan(err_check)), T);
        %     fprintf('  DEBUG: state_est(1,1)=%.4f, state_est(1,end)=%.4f\n', ...
        %         state_est(1,1), state_est(1,end));
        %     fprintf('  DEBUG: true_offset(1)=%.4f, true_offset(end)=%.4f\n', ...
        %         true_trajectories{i}.offset(1), true_trajectories{i}.offset(end));
        % 
        %     % Check dt_hours directly
        %     ts = cs.logger.timestamps(:)';
        %     is_v = ~isnat(ts);
        %     v_idx = find(is_v);
        %     dt_h = zeros(1, T);
        %     for kk = 2:numel(v_idx)
        %         dt_h(v_idx(kk)) = hours(ts(v_idx(kk)) - ts(v_idx(kk-1)));
        %     end
        %     dt_h(v_idx(1)) = 1;
        %     dt_h = min(dt_h, 48);
        %     total_t = sum(dt_h(is_v));
        %     tw = sum(dt_h(is_v) .* abs(err_check(is_v))) / total_t;
        %     fprintf('  DEBUG: manual tw_mae=%.6f, total_time=%.1f\n', tw, total_t);
        %     fprintf('  DEBUG: any NaN in dt_h(is_v)? %d\n', any(isnan(dt_h(is_v))));
        %     fprintf('  DEBUG: any NaN in err_check(is_v)? %d\n', any(isnan(err_check(is_v))));
        % end

        tw_mae_offset(i) = kpis.tw_mae_offset;
        tw_mae_U(i) = kpis.tw_mae_U;
        conv_days_offset(i) = kpis.convergence_days_offset;
        conv_days_U(i) = kpis.convergence_days_U;
        ss_mae_offset(i) = kpis.steady_state_mae_offset;
        ss_mae_U(i) = kpis.steady_state_mae_U;
        stability_offset(i) = kpis.stability_index_offset;
        stability_U(i) = kpis.stability_index_U;
        max_exc_offset(i) = kpis.max_excursion_offset;
        max_exc_U(i) = kpis.max_excursion_U;
        final_err_offset(i) = kpis.final_error_offset;
        final_err_U(i) = kpis.final_error_U;

        true_offset_initial(i) = true_trajectories{i}.offset(1);
        true_offset_final(i) = true_trajectories{i}.offset(end);
        true_offset_drift(i) = true_offset_final(i) - true_offset_initial(i);
        true_U_initial(i) = true_trajectories{i}.U(1);
        true_U_final(i) = true_trajectories{i}.U(end);

        initial_offsets(i) = kpis.initial_state.offset;
        initial_U(i) = kpis.initial_state.U;
        initial_std_offset(i) = kpis.initial_state.std_offset;
        initial_std_U(i) = kpis.initial_state.std_U;
        final_offsets(i) = kpis.final_state.offset;
        final_U_vals(i) = kpis.final_state.U;
        final_std_offset(i) = kpis.final_state.std_offset;
        final_std_U(i) = kpis.final_state.std_U;

        if ~isempty(kpis.end_of_season_errors)
            eos = kpis.end_of_season_errors;
            eos.house_id = repmat(house_ids(i), height(eos), 1);
            eos.csac_id = repmat(csac_id, height(eos), 1);
            all_eos = [all_eos; eos]; %#ok<AGROW>
        end
    end

    %% Build summary table
    summary = table(house_ids, ...
        true_offset_initial, true_offset_final, true_offset_drift, ...
        true_U_initial, true_U_final, ...
        tw_mae_offset, tw_mae_U, ...
        conv_days_offset, conv_days_U, ...
        ss_mae_offset, ss_mae_U, ...
        stability_offset, stability_U, ...
        max_exc_offset, max_exc_U, ...
        final_err_offset, final_err_U);
    summary.csac_id = repmat(csac_id, num_houses, 1);

    %% Save CSVs
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    writetable(summary, fullfile(output_folder, sprintf('kpis_csac_%d.csv', csac_id)));

    if ~isempty(all_eos)
        writetable(all_eos, fullfile(output_folder, sprintf('end_of_season_errors_csac_%d.csv', csac_id)));
    end

    %% Print summary
    fprintf('\n=== CSAC %d KPI Summary ===\n', csac_id);
    fprintf('  TW-MAE offset:        mean=%.3f °C,    median=%.3f °C\n', ...
        mean(tw_mae_offset, 'omitnan'), median(tw_mae_offset, 'omitnan'));
    fprintf('  TW-MAE U:             mean=%.4f,       median=%.4f W/m/K\n', ...
        mean(tw_mae_U, 'omitnan'), median(tw_mae_U, 'omitnan'));
    fprintf('  Convergence offset:   mean=%.0f days    (%d/%d converged)\n', ...
        mean(conv_days_offset, 'omitnan'), sum(isfinite(conv_days_offset)), num_houses);
    fprintf('  Convergence U:        mean=%.0f days    (%d/%d converged)\n', ...
        mean(conv_days_U, 'omitnan'), sum(isfinite(conv_days_U)), num_houses);
    fprintf('  SS-MAE offset:        mean=%.3f °C\n', mean(ss_mae_offset, 'omitnan'));
    fprintf('  SS-MAE U:             mean=%.4f W/m/K\n', mean(ss_mae_U, 'omitnan'));
    fprintf('  Stability offset:     mean=%.3f °C\n', mean(stability_offset, 'omitnan'));
    fprintf('  Stability U:          mean=%.4f W/m/K\n', mean(stability_U, 'omitnan'));
    fprintf('  Max excursion offset: mean=%.3f °C\n', mean(max_exc_offset, 'omitnan'));
    fprintf('  Max excursion U:      mean=%.4f W/m/K\n', mean(max_exc_U, 'omitnan'));
    fprintf('  Final error offset:   mean=%.3f °C\n', mean(abs(final_err_offset), 'omitnan'));
    fprintf('  Final error U:        mean=%.4f W/m/K\n', mean(abs(final_err_U), 'omitnan'));
    fprintf('  True offset drift:    mean=%.3f °C,    max=%.3f °C\n', ...
        mean(true_offset_drift, 'omitnan'), max(abs(true_offset_drift)));
    if ~isempty(all_eos)
        fprintf('  End-of-season errors (all houses):\n');
        for s = unique(all_eos.season)'
            s_data = all_eos(all_eos.season == s, :);
            fprintf('    Season %d: offset MAE=%.3f °C, U MAE=%.4f W/m/K\n', ...
                s, mean(abs(s_data.offset_error_C)), mean(abs(s_data.U_error_WmK)));
        end
    end
    fprintf('===========================\n\n');

    %% Generate visualization
    plot_kpi_bar_charts(house_ids, csac_id, ...
        true_offset_initial, true_offset_final, ...
        true_U_initial, true_U_final, ...
        initial_offsets, initial_U, initial_std_offset, initial_std_U, ...
        final_offsets, final_U_vals, final_std_offset, final_std_U, ...
        output_folder);
end