function save_diagnostic_summary(logger, ground_truth, csac_id, output_folder, filter_name)
%SAVE_DIAGNOSTIC_SUMMARY Exports a compact summary for external analysis.
%   Produces a small CSV with per-house final estimates, errors, and
%   key statistics about jumps and convergence.

    num_houses = size(logger.state_estimates, 2);
    num_timesteps = size(logger.state_estimates, 3);
    house_ids = logger.house_ids;
    
    % --- Per-house summary ---
    summary = table();
    for i = 1:num_houses
        offset_trace = squeeze(logger.state_estimates(1, i, :));
        U_trace = squeeze(logger.state_estimates(2, i, :));
        
        % Remove NaN/zero-padding from unused timesteps
        valid = offset_trace ~= 0 | U_trace ~= 0;
        if ~any(valid)
            continue
        end
        offset_trace = offset_trace(valid);
        U_trace = U_trace(valid);
        
        % Find the ground truth for this house
        gt_row = ground_truth(ground_truth.house_id == house_ids(i), :);
        if isempty(gt_row)
            true_offset = NaN;
            true_U = NaN;
        else
            true_offset = gt_row.true_offset;
            true_U = gt_row.true_U;
        end
        
        % Final estimates
        final_offset = offset_trace(end);
        final_U = U_trace(end);
        
        % Errors
        offset_error = final_offset - true_offset;
        U_error = final_U - true_U;
        
        % Detect largest single-step jumps
        offset_diffs = diff(offset_trace);
        U_diffs = diff(U_trace);
        [max_offset_jump, jump_idx_offset] = max(abs(offset_diffs));
        [max_U_jump, jump_idx_U] = max(abs(U_diffs));
        
        % Convergence: std of last 10% of trace
        tail_start = max(1, round(0.9 * length(offset_trace)));
        offset_tail_std = std(offset_trace(tail_start:end));
        U_tail_std = std(U_trace(tail_start:end));
        
        % RMSE over last 25% (steady-state accuracy)
        ss_start = max(1, round(0.75 * length(offset_trace)));
        offset_rmse = sqrt(mean((offset_trace(ss_start:end) - true_offset).^2));
        U_rmse = sqrt(mean((U_trace(ss_start:end) - true_U).^2));
        
        row = table(house_ids(i), final_offset, final_U, true_offset, true_U, ...
            offset_error, U_error, ...
            max_offset_jump, jump_idx_offset, ...
            max_U_jump, jump_idx_U, ...
            offset_tail_std, U_tail_std, ...
            offset_rmse, U_rmse, ...
            'VariableNames', {'house_id', 'final_offset', 'final_U', ...
            'true_offset', 'true_U', 'offset_error', 'U_error', ...
            'max_offset_jump', 'jump_idx_offset', ...
            'max_U_jump', 'jump_idx_U', ...
            'offset_tail_std', 'U_tail_std', ...
            'offset_rmse', 'U_rmse'});
        
        summary = [summary; row]; %#ok<AGROW>
    end
    
    % Save
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    filename = fullfile(output_folder, sprintf('%s_csac_%d_summary.csv', filter_name, csac_id));
    writetable(summary, filename);
    fprintf('Saved diagnostic summary: %s\n', filename);
end