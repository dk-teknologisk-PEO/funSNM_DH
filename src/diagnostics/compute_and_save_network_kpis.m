function summary = compute_and_save_network_kpis(cs, csac_id, output_folder, kpi_config)
%COMPUTE_AND_SAVE_NETWORK_KPIS Computes KPIs for all houses in a CSAC and saves results.
%
%   Builds constant true-value trajectories from the ground truth table.
%   When drift is introduced later, pass pre-built trajectories instead.
%
%   Args:
%       cs (struct): CSAC state struct after processing.
%       csac_id (int): CSAC identifier.
%       output_folder (char): Output directory for CSV files.
%       kpi_config (struct): KPI parameters (tolerances, hold days).
%
%   Returns:
%       summary (table): Summary table with one row per house.

    num_houses = cs.num_houses;
    T = length(cs.logger.timestamps);

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
    true_offset_val = nan(num_houses, 1);
    true_U_val = nan(num_houses, 1);

    all_eos = table();

    for i = 1:num_houses
        state_est = squeeze(cs.logger.state_estimates(:, i, :)); % 2xT

        % Build constant trajectories from ground truth
        true_offset_traj = repmat(cs.ground_truth.true_offset(i), 1, T);
        true_U_traj = repmat(cs.ground_truth.true_U(i), 1, T);

        kpis = compute_house_kpis(state_est, cs.logger.timestamps, ...
            true_offset_traj, true_U_traj, kpi_config);

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
        true_offset_val(i) = cs.ground_truth.true_offset(i);
        true_U_val(i) = cs.ground_truth.true_U(i);

        % Append end-of-season errors with house_id
        if ~isempty(kpis.end_of_season_errors)
            eos = kpis.end_of_season_errors;
            eos.house_id = repmat(house_ids(i), height(eos), 1);
            eos.csac_id = repmat(csac_id, height(eos), 1);
            all_eos = [all_eos; eos]; %#ok<AGROW>
        end
    end

    summary = table(house_ids, true_offset_val, true_U_val, ...
        tw_mae_offset, tw_mae_U, ...
        conv_days_offset, conv_days_U, ...
        ss_mae_offset, ss_mae_U, ...
        stability_offset, stability_U, ...
        max_exc_offset, max_exc_U);
    summary.csac_id = repmat(csac_id, num_houses, 1);

    % Save
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    writetable(summary, fullfile(output_folder, sprintf('kpis_csac_%d.csv', csac_id)));

    if ~isempty(all_eos)
        writetable(all_eos, fullfile(output_folder, sprintf('end_of_season_errors_csac_%d.csv', csac_id)));
    end

    % Print summary
    fprintf('\n=== CSAC %d KPI Summary ===\n', csac_id);
    fprintf('  TW-MAE offset:       mean=%.3f °C,   median=%.3f °C\n', ...
        mean(tw_mae_offset, 'omitnan'), median(tw_mae_offset, 'omitnan'));
    fprintf('  TW-MAE U:            mean=%.4f,      median=%.4f W/m/K\n', ...
        mean(tw_mae_U, 'omitnan'), median(tw_mae_U, 'omitnan'));
    fprintf('  Convergence offset:  mean=%.0f days   (%d/%d converged)\n', ...
        mean(conv_days_offset, 'omitnan'), sum(isfinite(conv_days_offset)), num_houses);
    fprintf('  Convergence U:       mean=%.0f days   (%d/%d converged)\n', ...
        mean(conv_days_U, 'omitnan'), sum(isfinite(conv_days_U)), num_houses);
    fprintf('  SS-MAE offset:       mean=%.3f °C\n', mean(ss_mae_offset, 'omitnan'));
    fprintf('  SS-MAE U:            mean=%.4f W/m/K\n', mean(ss_mae_U, 'omitnan'));
    fprintf('  Stability offset:    mean=%.3f °C\n', mean(stability_offset, 'omitnan'));
    fprintf('  Stability U:         mean=%.4f W/m/K\n', mean(stability_U, 'omitnan'));
    fprintf('  Max excursion offset: mean=%.3f °C\n', mean(max_exc_offset, 'omitnan'));
    fprintf('  Max excursion U:      mean=%.4f W/m/K\n', mean(max_exc_U, 'omitnan'));
    if ~isempty(all_eos)
        fprintf('  End-of-season errors (all houses):\n');
        for s = unique(all_eos.season)'
            s_data = all_eos(all_eos.season == s, :);
            fprintf('    Season %d: offset mean=%.3f °C, U mean=%.4f W/m/K\n', ...
                s, mean(abs(s_data.offset_error_C)), mean(abs(s_data.U_error_WmK)));
        end
    end
    fprintf('===========================\n\n');
end