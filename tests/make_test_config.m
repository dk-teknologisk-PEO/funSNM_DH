function config = make_test_config()
%MAKE_TEST_CONFIG Creates a minimal config struct for testing.

    config.project.cutoff.alpha_min = 0.3;
    config.project.cutoff.flow_cutoff = 10;
    config.project.cutoff.delta_T_gate_threshold = 5;
    config.project.cutoff.offset_min = -2.0;
    config.project.cutoff.offset_max = 2.0;
    config.project.cutoff.U_min = 0.05;
    config.project.cutoff.U_max = 0.25;

    config.project.initialization.min_active_houses = 3;
    config.project.initialization.max_innovation_C = 2.0;
    config.project.initialization.max_nis = 5.0;
    config.project.initialization.max_delta_T_change_rate = 5.0;
    config.project.initialization.innovation_gate_N_sigma = 3.0;

    config.project.initialization.ukf.measurement_noise = 0.1;
    config.project.initialization.ukf.process_noise_offset = 0.001;
    config.project.initialization.ukf.process_noise_U = 0.0001;
    config.project.initialization.ukf.state_uncertainty_offset = 2.0;
    config.project.initialization.ukf.state_uncertainty_U = 0.05;
    config.project.initialization.ukf.alpha_forget = 1.002;

    config.project.master_offset.gamma = 0.1;
    config.project.master_offset.deadzone = 0.01;

    config.project.debug.disable_stability_gate = false;
    config.project.debug.disable_mean_centering = false;
    config.project.debug.print_update_result = false;

    config.project.heating_season_gate.T_air_max_threshold = 8;
    config.project.heating_season_gate.lookback_days = 10;
    config.project.heating_season_gate.min_cooldown_days = 90;
    config.project.heating_season_gate.min_season_days_for_cooldown = 45;
    config.project.heating_season_gate.P_restart_fraction = 0.5;
    config.project.heating_season_gate.gap_P_growth_per_day_offset = 0.01;
    config.project.heating_season_gate.gap_P_growth_per_day_U = 0.001;

    config.project.simulation.offset_mean = 0.5;
    config.project.simulation.offset_std = 0.3;
    config.project.simulation.offset_min = 0.0;
    config.project.simulation.offset_max = 2.0;
    config.project.simulation.U_init_mean = 0.12;
    config.project.simulation.U_init_std = 0.01;
end