function plot_csac_U_diagnostics(all_cs, csac_ids, U_csac, network_id, output_folder)
%PLOT_CSAC_U_DIAGNOSTICS Plots offset and U_service vs position for all CSACs.

    num_csacs = numel(csac_ids);

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    fig = figure('Visible', 'off', 'Position', [50, 50, 600*num_csacs, 500]);

    for c = 1:num_csacs
        cs = all_cs{c};
        csac_id = csac_ids(c);
        num_houses = cs.num_houses;

        offsets = zeros(num_houses, 1);
        U_services = zeros(num_houses, 1);
        positions = zeros(num_houses, 1);
        std_offsets = zeros(num_houses, 1);
        std_U = zeros(num_houses, 1);

        for i = 1:num_houses
            offsets(i) = cs.ukf_states{i}.x(1);
            U_services(i) = cs.ukf_states{i}.x(2);
            std_offsets(i) = sqrt(abs(cs.ukf_states{i}.P(1,1)));
            std_U(i) = sqrt(abs(cs.ukf_states{i}.P(2,2)));
            house_id = cs.house_ids(i);
            house_rows = cs.meter_data(cs.meter_data.house_id == house_id, :);
            if ~isempty(house_rows)
                positions(i) = house_rows.x_pos_m(1);
            end
        end

        true_offsets = cs.ground_truth.true_offset;
        true_U = cs.ground_truth.true_U;

        % De-mean offsets for trend visualization
        offsets_dm = offsets - mean(offsets);
        true_offsets_dm = true_offsets - mean(true_offsets);

        % Fit trend
        p_off = polyfit(positions, offsets_dm, 1);
        x_fit = linspace(min(positions), max(positions), 100);

        %% Offset subplot
        subplot(2, num_csacs, c);
        hold on;
        errorbar(positions, offsets_dm, 2*std_offsets, 'bo', ...
            'MarkerFaceColor', [0.2 0.5 0.8], 'MarkerSize', 6, ...
            'LineWidth', 1, 'DisplayName', 'Estimated (de-meaned)');
        plot(positions, true_offsets_dm, 'r*', 'MarkerSize', 8, ...
            'LineWidth', 1.5, 'DisplayName', 'True (de-meaned)');
        plot(x_fit, polyval(p_off, x_fit), 'b--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Slope=%.4f', p_off(1)));
        yline(0, 'k:', 'HandleVisibility', 'off');
        hold off;
        xlabel('Position [m]');
        ylabel('Offset (de-meaned) [°C]');
        title(sprintf('CSAC %d — Offset vs Position', csac_id));
        legend('Location', 'best', 'FontSize', 7);
        grid on;

        %% U_service subplot
        subplot(2, num_csacs, num_csacs + c);
        hold on;
        errorbar(positions, U_services, 2*std_U, 'bo', ...
            'MarkerFaceColor', [0.8 0.5 0.2], 'MarkerSize', 6, ...
            'LineWidth', 1, 'DisplayName', 'Estimated');
        plot(positions, true_U, 'r*', 'MarkerSize', 8, ...
            'LineWidth', 1.5, 'DisplayName', 'True');
        hold off;
        xlabel('Position [m]');
        ylabel('U_{service} [W/m/K]');
        title(sprintf('CSAC %d — U_{service} vs Position', csac_id));
        legend('Location', 'best', 'FontSize', 7);
        grid on;
    end

    sgtitle(sprintf('Network %d — U_{csac}=%.4f (true=%.4f)', ...
        network_id, U_csac, all_cs{1}.U_csac_true));

    saveas(fig, fullfile(output_folder, sprintf('csac_U_diag_network_%d.png', network_id)));
    close(fig);
end