function [config, params] = load_config(config_file)
%LOAD_CONFIG Load config.json and extract all parameters into a flat struct.

    config = jsondecode(fileread(config_file));
    
    % Debug settings
    params.debug_disable_stability_gate = config.project.debug.disable_stability_gate;
    params.debug_disable_mean_centering = config.project.debug.disable_mean_centering;
    params.debug_print_update_result = config.project.debug.print_update_result;
    
    % Heating season gate
    params.season_gate_T_air_max_threshold = config.project.heating_season_gate.T_air_max_threshold;
    params.season_gate_lookback_days = config.project.heating_season_gate.lookback_days;
    params.season_gate_P_restart_fraction = config.project.heating_season_gate.P_restart_fraction;
    params.season_gate_min_cooldown_days = config.project.heating_season_gate.min_cooldown_days;
    params.season_gate_min_season_for_cooldown = config.project.heating_season_gate.min_season_days_for_cooldown;
    
    % Master offset
    params.master_offset_gamma = config.project.master_offset.gamma;
    params.master_offset_deadzone = config.project.master_offset.deadzone;
    
    % Innovation / NIS gates
    params.hard_innovation_gate_C = config.project.initialization.max_innovation_C;
    params.max_nis = config.project.initialization.max_nis;
    
    % UKF base matrices
    params.R_base = config.project.initialization.ukf.measurement_noise^2;
    params.Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, ...
                          (config.project.initialization.ukf.process_noise_U)^2]);
    params.P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, ...
                          (config.project.initialization.ukf.state_uncertainty_U)^2]);
    params.innovation_gate_initial = config.project.initialization.ukf.state_uncertainty_offset * ...
                                     config.project.initialization.innovation_gate_N_sigma;
    params.P_max_restart = diag([config.project.heating_season_gate.P_max_offset_restart^2, ...
                                 config.project.heating_season_gate.P_max_U_restart^2]);
    
    % PF
    params.num_particles = config.project.initialization.pf.num_particles;
    
    % Gating
    params.absolute_flow_floor_kg_h = config.project.cutoff.flow_cutoff;
    params.delta_T_gate_threshold = config.project.cutoff.delta_T_gate_threshold;
    params.min_active_houses = config.project.initialization.min_active_houses;
    params.alpha_min = config.project.cutoff.alpha_min;
    params.max_delta_T_change = config.project.initialization.max_delta_T_change_rate;
    
    % Mean centering
    params.error_meaning = config.project.initialization.modify_mean_offset;
    
    % Pass config through for functions that need the full struct
    params.config = config;
    
    fprintf('Hard innovation gate set to %.2f °C\n', params.hard_innovation_gate_C);
end