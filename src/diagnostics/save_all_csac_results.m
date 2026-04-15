function save_all_csac_results(csac_state, csac, network, output_folder_ukf, output_folder_pf)
%SAVE_ALL_CSAC_RESULTS Save all diagnostic outputs for one CSAC.

    gt = csac_state.ground_truth;
    
    plot_diagnostics(csac_state.logger_ukf, gt, csac, network, output_folder_ukf);
    plot_diagnostics(csac_state.logger_pf, gt, csac, network, output_folder_pf);
    
    save_logger_to_csv(csac_state.logger_ukf, output_folder_ukf, strcat('ukf_csac_', string(csac)));
    save_logger_to_csv(csac_state.logger_pf, output_folder_pf, strcat('pf_csac_', string(csac)));
    
    save_diagnostic_summary(csac_state.logger_ukf, gt, csac, output_folder_ukf, 'ukf');
    save_diagnostic_summary(csac_state.logger_pf, gt, csac, output_folder_pf, 'pf');
    
    save_diagnostic_summary_detailed(csac_state.logger_ukf, gt, csac, output_folder_ukf, 'ukf');
    save_diagnostic_summary_detailed(csac_state.logger_pf, gt, csac, output_folder_pf, 'pf');
    
    save_daily_diagnostics(csac_state.logger_ukf, gt, csac, output_folder_ukf, 'ukf');
    save_daily_diagnostics(csac_state.logger_pf, gt, csac, output_folder_pf, 'pf');
end