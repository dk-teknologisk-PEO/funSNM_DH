function plot_kpi_bar_charts(house_ids, csac_id, ...
    true_offset_initial, true_offset_final, ...
    true_U_initial, true_U_final, ...
    init_offset, init_U, init_std_offset, init_std_U, ...
    final_offset, final_U, final_std_offset, final_std_U, ...
    output_folder)

    % Guard against complex values from negative covariance
    init_offset = real(init_offset);
    init_U = real(init_U);
    init_std_offset = real(init_std_offset);
    init_std_U = real(init_std_U);
    final_offset = real(final_offset);
    final_U = real(final_U);
    final_std_offset = real(final_std_offset);
    final_std_U = real(final_std_U);
    true_offset_initial = real(true_offset_initial);
    true_offset_final = real(true_offset_final);
    true_U_initial = real(true_U_initial);
    true_U_final = real(true_U_final);
%PLOT_KPI_BAR_CHARTS Generates bar charts comparing initial and final estimates.
%
%   Offset plot:
%     - initial estimate ± 2σ
%     - final estimate ± 2σ
%     - true initial value
%     - true final value
%
%   U-value plot:
%     - initial estimate ± 2σ
%     - final estimate ± 2σ
%     - true initial/final values (typically identical)

    n = numel(house_ids);
    x = 1:n;
    house_labels = arrayfun(@(id) sprintf('H%d', id), house_ids, 'UniformOutput', false);

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    %% ============================================================
    %% Figure 1: Offset Estimates
    %% ============================================================
    fig1 = figure('Visible', 'off', 'Position', [100, 100, 1000, 500]);

    bar_width = 0.35;

    bar(x - bar_width/2, init_offset, bar_width, ...
        'FaceColor', [0.7 0.7 0.9], 'DisplayName', 'Initial estimate');
    hold on;
    bar(x + bar_width/2, final_offset, bar_width, ...
        'FaceColor', [0.2 0.5 0.8], 'DisplayName', 'Final estimate');

    errorbar(x - bar_width/2, init_offset, 2*init_std_offset, 'k.', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    errorbar(x + bar_width/2, final_offset, 2*final_std_offset, 'k.', ...
        'LineWidth', 1, 'HandleVisibility', 'off');

    plot(x, true_offset_initial, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5, ...
        'DisplayName', 'True initial');
    plot(x, true_offset_final, 'r*', 'MarkerSize', 9, 'LineWidth', 1.5, ...
        'DisplayName', 'True final');

    set(gca, 'XTick', x, 'XTickLabel', house_labels);
    xlabel('House');
    ylabel('Offset [°C]');
    title(sprintf('CSAC %d — Offset: Initial vs Final Estimate (95%% CI)', csac_id));
    legend('Location', 'best');
    grid on;
    hold off;

    saveas(fig1, fullfile(output_folder, sprintf('kpi_offset_bar_csac_%d.png', csac_id)));
    close(fig1);

    %% ============================================================
    %% Figure 2: U-Value Estimates
    %% ============================================================
    fig2 = figure('Visible', 'off', 'Position', [100, 100, 1000, 500]);

    bar(x - bar_width/2, init_U, bar_width, ...
        'FaceColor', [0.9 0.8 0.7], 'DisplayName', 'Initial estimate');
    hold on;
    bar(x + bar_width/2, final_U, bar_width, ...
        'FaceColor', [0.8 0.5 0.2], 'DisplayName', 'Final estimate');

    errorbar(x - bar_width/2, init_U, 2*init_std_U, 'k.', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    errorbar(x + bar_width/2, final_U, 2*final_std_U, 'k.', ...
        'LineWidth', 1, 'HandleVisibility', 'off');

    plot(x, true_U_initial, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5, ...
        'DisplayName', 'True initial');
    plot(x, true_U_final, 'r*', 'MarkerSize', 9, 'LineWidth', 1.5, ...
        'DisplayName', 'True final');

    set(gca, 'XTick', x, 'XTickLabel', house_labels);
    xlabel('House');
    ylabel('U-value [W/m/K]');
    title(sprintf('CSAC %d — U-value: Initial vs Final Estimate (95%% CI)', csac_id));
    legend('Location', 'best');
    grid on;
    hold off;

    saveas(fig2, fullfile(output_folder, sprintf('kpi_U_bar_csac_%d.png', csac_id)));
    close(fig2);
end