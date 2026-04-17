function csac = initialize_csac_state(topology, topology_csac, houses_csac, meter_data, num_timestamps, Q_base, R_base, P_base, innovation_gate_initial)
%INITIALIZE_CSAC_STATE Creates all per-CSAC state needed for the main loop.
%
%   Packages filter states, snapshots, loggers, gate counters, and
%   precomputed data tables into a single struct.
%
%   Args:
%       topology (struct): Full network topology.
%       topology_csac (struct array): Topology entries for houses in this CSAC.
%       houses_csac (logical): Mask of houses belonging to this CSAC.
%       meter_data (table): Full meter data table for the network.
%       num_timestamps (int): Total number of timesteps for logger allocation.
%       Q_base (2x2 matrix): Base process noise covariance.
%       R_base (scalar): Base measurement noise variance.
%       P_base (2x2 matrix): Initial state covariance.
%       innovation_gate_initial (scalar): Initial innovation gate value.
%
%   Returns:
%       csac (struct): Complete CSAC state struct with fields:
%           .num_houses (int): Number of houses in the CSAC.
%           .house_ids (1xN int): Sorted house IDs.
%           .U_csac (scalar): CSAC pipe U-value from topology.
%           .ground_truth (table): True offset and U-value per house.
%           .meter_data (table): Meter data joined with geometry, for this CSAC.
%           .ukf_states (cell): UKF state structs per house.
%           .ukf_states_snapshot (cell): Snapshot of UKF states per house.
%           .ukf_innovation_gate (1xN double): Innovation gate per house.
%           .ukf_offsets (Nx1 double): Current offset estimates (NaN until first update).
%           .last_valid_T_main_C (Nx1 double): Last valid T_main per house.
%           .last_update_timestamp (Nx1 datetime): Last update time per house.
%           .gate_reject_count (int): Total rejected updates.
%           .gate_accept_count (int): Total accepted updates.
%           .season_state (struct): Heating season gate state.
%           .logger (struct): UKF logger struct.

    num_houses = sum(houses_csac);
    house_ids = sort([topology.houses(houses_csac).id]);

    %% Ground truth table
    true_offset = zeros(num_houses, 1);
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

    %% Initialize UKF states
    ukf_states = cell(1, num_houses);
    ukf_states_snapshot = cell(1, num_houses);
    ukf_innovation_gate = nan(1, num_houses);

    for i = 1:num_houses
        x_init = [randn() * 0.3; 0.12 + randn() * 0.03];
        ukf_states{i} = initialize_kalman_filter(x_init, Q_base, R_base, P_base);
        ukf_innovation_gate(i) = innovation_gate_initial;
        ukf_states_snapshot{i}.x = ukf_states{i}.x;
        ukf_states_snapshot{i}.P = ukf_states{i}.P;
    end

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
    csac.season_state = initialize_season_state();
    csac.logger = logger;
end