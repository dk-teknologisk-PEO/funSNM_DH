function plot_diagnostics(logger, ground_truth, csac_id, network_name, output_folder)
%PLOT_DIAGNOSTICS Generates and saves a comprehensive diagnostic figure.
%   Highlights regions where |error| > 2*sigma (95% CI doesn't cover truth).

    fig_name = strcat('Diagnostics_Network_CSAC_', string(network_name), '_', string(csac_id));
    fig = figure('Name', fig_name, 'Position', [50, 50, 1600, 900], 'Visible', 'off');

    num_houses = length(logger.house_ids);
    timestamps = logger.timestamps;
    colors = lines(num_houses);

    %% Plot 1: Offset Estimation Error with uncovered highlights
    subplot(3, 2, 1);
    hold on;
    for i = 1:num_houses
        err = squeeze(logger.state_estimates(1, i, :) - ground_truth.true_offset(i));
        std_val = squeeze(sqrt(abs(logger.covariance_posterior(1, i, :))));
        valid = ~isnat(timestamps') & ~isnan(err);

        % Base line
        plot(timestamps, err, '-', 'Color', colors(i,:), 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');

        % Highlight uncovered points
        uncovered = valid & (abs(err) > 2 * std_val);
        if any(uncovered)
            plot(timestamps(uncovered), err(uncovered), 'o', ...
                'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], ...
                'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end
    end
    % Legend entries (separate loop for clean legend)
    for i = 1:num_houses
        plot(NaT, NaN, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('House %d', logger.house_ids(i)));
    end
    plot(NaT, NaN, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0 0], ...
        'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.5, ...
        'DisplayName', '|err| > 2\sigma');
    yline(0, 'k--', 'HandleVisibility', 'off');
    hold off;
    title('Offset Estimation Error (red = outside 95% CI)');
    ylabel('[°C]'); grid on;
    legend('Location', 'best', 'FontSize', 7);

    %% Plot 2: U-Value Estimation Error with uncovered highlights
    subplot(3, 2, 2);
    hold on;
    for i = 1:num_houses
        err = squeeze(logger.state_estimates(2, i, :) - ground_truth.true_U(i));
        std_val = squeeze(sqrt(logger.covariance_posterior(2, i, :)));
        valid = ~isnat(timestamps') & ~isnan(err);

        plot(timestamps, err, '-', 'Color', colors(i,:), 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');

        uncovered = valid & (abs(err) > 2 * std_val);
        if any(uncovered)
            plot(timestamps(uncovered), err(uncovered), '.', ...
                'Color', [1 0 0], 'MarkerSize', 4, 'HandleVisibility', 'off');
        end
    end
    for i = 1:num_houses
        plot(NaT, NaN, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('House %d', logger.house_ids(i)));
    end
    plot(NaT, NaN, '.', 'Color', [1 0 0], 'MarkerSize', 8, ...
        'DisplayName', '|err| > 2\sigma');
    yline(0, 'k--', 'HandleVisibility', 'off');
    hold off;
    title('U-Value Estimation Error (red = outside 95% CI)');
    ylabel('[W/m/K]'); grid on;
    legend('Location', 'best', 'FontSize', 7);

    %% Plot 3: Innovation (Residual) with 3-sigma bounds
    subplot(3, 2, 3);
    hold on;
    i = 1;
    sigma_bounds = 3 * sqrt(squeeze(logger.innovation_variance(i, :)));
    plot(timestamps, squeeze(logger.innovations(i, :)), 'b.');
    plot(timestamps, sigma_bounds, 'r--', 'LineWidth', 1.5);
    plot(timestamps, -sigma_bounds, 'r--', 'LineWidth', 1.5);
    hold off;
    title(sprintf('Innovation (Residual) for House %d', logger.house_ids(i)));
    ylabel('Temp [°C]'); grid on; legend('Innovation', '3-Sigma Bound');

    %% Plot 4: Normalized Innovation Squared (NIS)
    subplot(3, 2, 4);
    hold on;
    nis_95_confidence = 3.84;
    for i = 1:num_houses
        plot(timestamps, squeeze(logger.nis(i, :)), '.', 'MarkerSize', 3);
    end
    yline(nis_95_confidence, 'r--', 'LineWidth', 1.5, 'Label', '95% Confidence Bound');
    hold off;
    title('Normalized Innovation Squared (NIS)');
    ylabel('NIS Value'); set(gca, 'YScale', 'log'); grid on;

    %% Plot 5: State Uncertainty - Offset
    subplot(3, 2, 5);
    hold on;
    for i = 1:num_houses
        stdev = sqrt(squeeze(logger.covariance_posterior(1, i, :)));
        plot(timestamps, stdev, '-', 'Color', colors(i,:), 'LineWidth', 1.0);
    end
    hold off;
    title('State Uncertainty (\sigma) - Offset');
    ylabel('Std. Dev. [°C]'); grid on;

    %% Plot 6: State Uncertainty - U-Value
    subplot(3, 2, 6);
    hold on;
    for i = 1:num_houses
        stdev = sqrt(squeeze(logger.covariance_posterior(2, i, :)));
        plot(timestamps, stdev, '-', 'Color', colors(i,:), 'LineWidth', 1.0);
    end
    hold off;
    title('State Uncertainty (\sigma) - U-Value');
    ylabel('Std. Dev. [W/m/K]'); grid on;

    %% Save
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    saveas(fig, fullfile(output_folder, strcat(fig_name, '.png')));
    close(fig);
    fprintf('Saved diagnostic plot to %s\n', fullfile(output_folder, strcat(fig_name, '.png')));
end