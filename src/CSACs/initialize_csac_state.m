function csac = initialize_csac_state(topology, topology_csac, houses_csac, ...
    meter_data, num_timestamps, Q_base, R_base, P_base, innovation_gate_initial, config, network_id)
%INITIALIZE_CSAC_STATE Creates all per-CSAC state needed for the main loop.

    num_houses = sum(houses_csac);
    house_ids = sort([topology.houses(houses_csac).id]);

    %% Read simulation config
    sim_cfg = config.project.simulation;

    %% Generate true offsets (reproducible per network + house)
    true_offset = zeros(num_houses, 1);
    for i = 1:num_houses
        rng(network_id * 10000 + house_ids(i), 'twister');
        raw_offset = sim_cfg.offset_mean + sim_cfg.offset_std * randn();
        true_offset(i) = max(sim_cfg.offset_min, min(sim_cfg.offset_max, raw_offset));
    end

    %% Ground truth table
    ground_truth = table(house_ids', true_offset, ...
        [topology_csac.service_pipe_insulation_W_m_K]', ...
        'VariableNames', {'house_id', 'true_offset', 'true_U'});

    %% Geometry table and meter data join
    csac_table = table(house_ids', ...
        [topology_csac.service_pipe_len_m]', ...
        [topology_csac.dist_on_cul_de_sac_m]', ...
        'VariableNames', {'house_id', 'length_service_m', 'x_pos_m'});

    meter_data_csac = meter_data(ismember(meter_data.house_id, house_ids), :);
    meter_data_csac = join(meter_data_csac, csac_table, 'Keys', 'house_id');

    %% Apply true offsets to meter data
    for i = 1:num_houses
        mask = meter_data_csac.house_id == house_ids(i);
        meter_data_csac.T_supply_C(mask) = meter_data_csac.T_supply_C(mask) - true_offset(i);
    end

    fprintf('  Applied simulated offsets to %d houses: [%s] °C\n', ...
        num_houses, strjoin(arrayfun(@(x) sprintf('%.2f', x), true_offset, 'UniformOutput', false), ', '));

    %% Initialize UKF states at prior (reproducible per network + house)
    ukf_states = cell(1, num_houses);
    ukf_states_snapshot = cell(1, num_houses);
    ukf_innovation_gate = nan(1, num_houses);

    for i = 1:num_houses
        rng(network_id * 10000 + house_ids(i) + 1e6, 'twister');
        x_init = [0; sim_cfg.U_init_mean + sim_cfg.U_init_std * randn()];
        ukf_states{i} = initialize_kalman_filter(x_init, Q_base, R_base, P_base);
        ukf_innovation_gate(i) = innovation_gate_initial;
        ukf_states_snapshot{i}.x = ukf_states{i}.x;
        ukf_states_snapshot{i}.P = ukf_states{i}.P;
    end

    %% Restore RNG
    rng('shuffle');

    %% Initialize logger and log initial states
    logger = initialize_logger(num_houses, num_timestamps, house_ids);
    for i = 1:num_houses
        logger.state_estimates(:, i, 1) = ukf_states{i}.x;
        logger.covariance_posterior(:, i, 1) = [ukf_states{i}.P(1,1); ukf_states{i}.P(2,2)];
    end

    %% Pack into struct
    csac.num_houses = num_houses;
    csac.house_ids = house_ids;
    csac.U_csac = topology.pipe_parameters.csac_pipe.insulation_W_m_K;
    csac.ground_truth = ground_truth;
    csac.meter_data = meter_data_csac;
    csac.ukf_states = ukf_states;
    csac.ukf_states_snapshot = ukf_states_snapshot;
    csac.ukf_innovation_gate = ukf_innovation_gate;
    csac.ukf_offsets = nan(num_houses, 1);
    csac.last_valid_T_main_C = nan(num_houses, 1);
    csac.last_update_timestamp = NaT(num_houses, 1);
    csac.gate_reject_count = 0;
    csac.gate_accept_count = 0;
    csac.consecutive_rejection_counters = zeros(num_houses, 1);
    csac.season_state = initialize_season_state();
    csac.logger = logger;
end