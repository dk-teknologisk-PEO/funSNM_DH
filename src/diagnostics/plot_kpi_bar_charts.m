function plot_kpi_bar_charts(house_ids, csac_id, ...
    true_offset, true_U, ...
    init_offset, init_U, init_std_offset, init_std_U, ...
    final_offset, final_U, final_std_offset, final_std_U, ...
    output_folder)
%PLOT_KPI_BAR_CHARTS Generates bar charts comparing initial and final estimates.
%
%   Creates two figures:
%   1. Offset estimates: initial vs final vs true, with 95% uncertainty bars
%   2. U-value estimates: initial vs final vs true, with 95% uncertainty bars

    n = numel(house_ids);
    x = 1:n;
    house_labels = arrayfun(@(id) sprintf('H%d', id), house_ids, 'UniformOutput', false);

    %% Figure 1: Offset Estimates
    fig1 = figure('Visible', 'off', 'Position', [100, 100, 900, 450]);

    bar_width = 0.35;
    bar(x - bar_width/2, init_offset, bar_width, 'FaceColor', [0.7 0.7 0.9], ...
        'DisplayName', 'Initial');
    hold on;
    bar(x + bar_width/2, final_offset, bar_width, 'FaceColor', [0.2 0.5 0.8], ...
        'DisplayName', 'Final');
    errorbar(x - bar_width/2, init_offset, 2*init_std_offset, 'k.', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    errorbar(x + bar_width/2, final_offset, 2*final_std_offset, 'k.', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    plot(x, true_offset, 'r*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'True');

    set(gca, 'XTick', x, 'XTickLabel', house_labels);
    xlabel('House');
    ylabel('Offset [°C]');
    title(sprintf('CSAC %d — Offset: Initial vs Final Estimate (95%% CI)', csac_id));
    legend('Location', 'best');
    grid on;
    hold off;

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    saveas(fig1, fullfile(output_folder, sprintf('kpi_offset_bar_csac_%d.png', csac_id)));

    %% Figure 2: U-Value Estimates
    fig2 = figure('Visible', 'off', 'Position', [100, 100, 900, 450]);

    bar(x - bar_width/2, init_U, bar_width, 'FaceColor', [0.9 0.8 0.7], ...
        'DisplayName', 'Initial');
    hold on;
    bar(x + bar_width/2, final_U, bar_width, 'FaceColor', [0.8 0.5 0.2], ...
        'DisplayName', 'Final');
    errorbar(x - bar_width/2, init_U, 2*init_std_U, 'k.', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    errorbar(x + bar_width/2, final_U, 2*final_std_U, 'k.', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    plot(x, true_U, 'r*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'True');

    set(gca, 'XTick', x, 'XTickLabel', house_labels);
    xlabel('House');
    ylabel('U-value [W/m/K]');
    title(sprintf('CSAC %d — U-value: Initial vs Final Estimate (95%% CI)', csac_id));
    legend('Location', 'best');
    grid on;
    hold off;

    saveas(fig2, fullfile(output_folder, sprintf('kpi_U_bar_csac_%d.png', csac_id)));
end