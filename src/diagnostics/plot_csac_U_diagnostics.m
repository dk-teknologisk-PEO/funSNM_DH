function plot_csac_U_diagnostics(cs, csac_id, network_id, output_folder)
%PLOT_CSAC_U_DIAGNOSTICS Plots diagnostics for CSAC pipe U-value estimation.
%
%   Creates three subplots:
%   1. Estimated offset vs position along CSAC (should show no trend if U_csac is correct)
%   2. Estimated U_service vs position (should show no trend)
%   3. U_csac evolution over time (convergence plot)

    num_houses = cs.num_houses;

    % Extract current estimates and positions
    offsets = zeros(num_houses, 1);
    U_service = zeros(num_houses, 1);
    positions = zeros(num_houses, 1);
    std_offsets = zeros(num_houses, 1);
    std_U = zeros(num_houses, 1);

    for i = 1:num_houses
        offsets(i) = cs.ukf_states{i}.x(1);
        U_service(i) = cs.ukf_states{i}.x(2);
        std_offsets(i) = sqrt(cs.ukf_states{i}.P(1,1));
        std_U(i) = sqrt(cs.ukf_states{i}.P(2,2));

        house_id = cs.house_ids(i);
        house_rows = cs.meter_data(cs.meter_data.house_id == house_id, :);
        if ~isempty(house_rows)
            positions(i) = house_rows.x_pos_m(1);
        end
    end

    % True values for comparison
    true_offsets = cs.ground_truth.true_offset;
    true_U_service = cs.ground_truth.true_U;

    house_labels = arrayfun(@(id) sprintf('H%d', id), cs.house_ids, 'UniformOutput', false);

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    fig = figure('Visible', 'off', 'Position', [100, 100, 1400, 900]);

    %% Subplot 1: Offset vs Position
    subplot(2, 2, 1);
    hold on;
    errorbar(positions, offsets, 2*std_offsets, 'bo', 'MarkerFaceColor', 'b', ...
        'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'Estimated');
    plot(positions, true_offsets, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'DisplayName', 'True');

    % Add linear fit line
    p = polyfit(positions, offsets, 1);
    x_fit = linspace(min(positions), max(positions), 100);
    plot(x_fit, polyval(p, x_fit), 'b--', 'LineWidth', 1, ...
        'DisplayName', sprintf('Fit: slope=%.4f °C/m', p(1)));

    % Label each point
    for i = 1:num_houses
        text(positions(i), offsets(i) + 0.05, house_labels{i}, ...
            'FontSize', 7, 'HorizontalAlignment', 'center');
    end

    xlabel('Position along CSAC [m]');
    ylabel('Offset [°C]');
    title('Estimated Offset vs Position');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    hold off;

    %% Subplot 2: U_service vs Position
    subplot(2, 2, 2);
    hold on;
    errorbar(positions, U_service, 2*std_U, 'bo', 'MarkerFaceColor', 'b', ...
        'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'Estimated');
    plot(positions, true_U_service, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'DisplayName', 'True');

    p_U = polyfit(positions, U_service, 1);
    plot(x_fit, polyval(p_U, x_fit), 'b--', 'LineWidth', 1, ...
        'DisplayName', sprintf('Fit: slope=%.6f /m', p_U(1)));

    for i = 1:num_houses
        text(positions(i), U_service(i) + 0.005, house_labels{i}, ...
            'FontSize', 7, 'HorizontalAlignment', 'center');
    end

    xlabel('Position along CSAC [m]');
    ylabel('U-value [W/m/K]');
    title('Estimated U_{service} vs Position');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    hold off;

    %% Subplot 3: U_csac evolution
    subplot(2, 2, 3);
    if isfield(cs, 'csac_U_history') && numel(cs.csac_U_history) > 1
        plot(1:numel(cs.csac_U_history), cs.csac_U_history, 'b-', 'LineWidth', 1.5, ...
            'DisplayName', 'Estimated');
        hold on;
        yline(cs.U_csac_true, 'r--', 'LineWidth', 1.5, 'DisplayName', 'True');
        xlabel('Update number');
        ylabel('U_{csac} [W/m/K]');
        title(sprintf('U_{csac} Convergence (final=%.4f, true=%.4f)', cs.U_csac, cs.U_csac_true));
        legend('Location', 'best');
        grid on;
        hold off;
    else
        text(0.5, 0.5, 'No U_{csac} history available', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        title('U_{csac} Convergence');
    end

    %% Subplot 4: Offset error vs Position
    subplot(2, 2, 4);
    hold on;
    offset_errors = offsets - true_offsets;
    bar(1:num_houses, offset_errors, 'FaceColor', [0.5 0.7 1.0]);

    set(gca, 'XTick', 1:num_houses, 'XTickLabel', house_labels);
    xlabel('House (ordered by position)');
    ylabel('Offset Error [°C]');
    title(sprintf('Offset Error by House (U_{csac}=%.4f, true=%.4f)', cs.U_csac, cs.U_csac_true));
    yline(0, 'k--');
    grid on;
    hold off;

    %% Add overall title
    sgtitle(sprintf('Network %d, CSAC %d — Pipe U-value Diagnostics', network_id, csac_id), ...
        'FontSize', 14);

    %% Save
    prefix = sprintf('network_%d_csac_%d', network_id, csac_id);
    save_figure(fig, fullfile(output_folder, sprintf('csac_U_diag_network_%d', network_id)));
    close(fig);
    fprintf('Saved CSAC U diagnostics to %s\n', ...
        fullfile(output_folder, sprintf('csac_U_diagnostics_%s.png', prefix)));
end