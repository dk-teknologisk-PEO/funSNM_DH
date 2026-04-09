% In a new file: plot_diagnostics.m

function plot_diagnostics(logger, ground_truth, csac_id, network_name, output_folder)
% PLOT_DIAGNOSTICS Generates and saves a comprehensive diagnostic figure.

    % Create a unique figure name
    fig_name = strcat('Diagnostics_Network_CSAC_', string(network_name),'_', string(csac_id));
    fig = figure('Name', fig_name, 'Position', [50, 50, 1600, 900], 'Visible', 'off');

    num_houses = length(logger.house_ids);
    timestamps = logger.timestamps;
    
    % Plot 1: Offset Error
    subplot(3, 2, 1);
    hold on;
    for i = 1:num_houses
        offset_error = logger.state_estimates(1, i, :) - ground_truth.true_offset(i);
        plot(timestamps, squeeze(offset_error), 'LineWidth', 1.5);
    end
    hold off;
    title('Offset Estimation Error (Estimated - True)');
    ylabel('[°C]'); grid on; legend(strcat('House ', string(logger.house_ids)));

    % Plot 2: U-Value Error
    subplot(3, 2, 2);
    hold on;
    for i = 1:num_houses
        U_error = logger.state_estimates(2, i, :) - ground_truth.true_U(i);
        plot(timestamps, squeeze(U_error), 'LineWidth', 1.5);
    end
    hold off;
    title('U-Value Estimation Error (Estimated - True)');
    ylabel('[W/m/K]'); grid on; legend(strcat('House ', string(logger.house_ids)));

    % Plot 3: Innovation (Residual) with 3-sigma bounds
    subplot(3, 2, 3);
    hold on;
    i = 1; % Plot for the first house as an example
    sigma_bounds = 3 * sqrt(squeeze(logger.innovation_variance(i, :)));
    plot(timestamps, squeeze(logger.innovations(i, :)), 'b.');
    plot(timestamps, sigma_bounds, 'r--', 'LineWidth', 1.5);
    plot(timestamps, -sigma_bounds, 'r--', 'LineWidth', 1.5);
    hold off;
    title(sprintf('Innovation (Residual) for House %d', logger.house_ids(i)));
    ylabel('Temp [°C]'); grid on; legend('Innovation', '3-Sigma Bound');

    % Plot 4: Normalized Innovation Squared (NIS)
    subplot(3, 2, 4);
    hold on;
    % Chi-squared 95% confidence interval for 1 DoF is ~3.84
    nis_95_confidence = 3.84;
    for i = 1:num_houses
        plot(timestamps, squeeze(logger.nis(i, :)), '.');
    end
    yline(nis_95_confidence, 'r--', 'LineWidth', 1.5, 'Label', '95% Confidence Bound');
    hold off;
    title('Normalized Innovation Squared (NIS)');
    ylabel('NIS Value'); set(gca, 'YScale', 'log'); grid on;

    % Plot 5: State Covariance (Uncertainty) - Offset
    subplot(3, 2, 5);
    hold on;
    for i = 1:num_houses
        stdev = sqrt(squeeze(logger.covariance_posterior(1, i, :)));
        plot(timestamps, stdev, 'LineWidth', 1.5);
    end
    hold off;
    title('State Uncertainty (sqrt(P)) - Offset');
    ylabel('Std. Dev. [°C]'); grid on;

    % Plot 6: State Covariance (Uncertainty) - U-Value
    subplot(3, 2, 6);
    hold on;
    for i = 1:num_houses
        stdev = sqrt(squeeze(logger.covariance_posterior(2, i, :)));
        plot(timestamps, stdev, 'LineWidth', 1.5);
    end
    hold off;
    title('State Uncertainty (sqrt(P)) - U-Value');
    ylabel('Std. Dev. [W/m/K]'); grid on;

    % Save the figure
    if ~exist(output_folder, 'dir')
       mkdir(output_folder)
    end
    saveas(fig, fullfile(output_folder, strcat(fig_name, '.png')), 'Resolution', 100);
    close(fig);
    fprintf('Saved diagnostic plot to %s\n',fullfile(output_folder, strcat(fig_name, '.png')));
end