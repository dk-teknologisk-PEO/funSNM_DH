function all_kpi_summaries = post_process_all_csacs(all_cs, all_true_traj, csac_ids, ...
    network_id, output_folder, kpi_config)
%POST_PROCESS_ALL_CSACS Runs post-processing for all CSACs after the main loop.
%
%   Args:
%       all_cs (cell array): CSAC state structs after processing.
%       all_true_traj (cell array): True trajectory cells per CSAC.
%       csac_ids (int array): CSAC identifiers.
%       network_id (int): Network identifier.
%       output_folder (char): Output directory.
%       kpi_config (struct): KPI configuration.
%
%   Returns:
%       all_kpi_summaries (table): Aggregated KPI summary across all CSACs.

    all_kpi_summaries = table();

    for c = 1:numel(csac_ids)
        csac_id = csac_ids(c);
        cs = all_cs{c};
        true_traj = all_true_traj{c};

        % Print statistics
        print_csac_summary(csac_id, cs.gate_accept_count, cs.gate_reject_count, cs.season_state);

        % Compute and save KPIs
        kpi_summary = compute_and_save_network_kpis(cs, csac_id, output_folder, ...
            kpi_config, true_traj);
        all_kpi_summaries = [all_kpi_summaries; kpi_summary]; %#ok<AGROW>

        % Plot diagnostics
        plot_diagnostics(cs.logger, cs.ground_truth, csac_id, network_id, output_folder);
        save_logger_to_csv(cs.logger, output_folder, strcat('ukf_csac_', string(csac_id)));
        save_diagnostic_summary(cs.logger, cs.ground_truth, csac_id, output_folder, 'ukf');
        save_diagnostic_summary_detailed(cs.logger, cs.ground_truth, csac_id, output_folder, 'ukf');
        save_daily_diagnostics(cs.logger, cs.ground_truth, csac_id, output_folder, 'ukf');

        % Close diagnostic figures but keep waitbar
        figs = findall(0, 'Type', 'figure');
        for fig = figs'
            if ~strcmp(get(fig, 'Tag'), 'TMWWaitbar')
                close(fig);
            end
        end
    end
end