function csac_state = initialize_csac(topology, meter_data, csac, params, timestamps)
%INITIALIZE_CSAC Set up all state variables for one cul-de-sac.

    % Extract topology
    U_csac = topology.pipe_parameters.csac_pipe.insulation_W_m_K;
    houses_csac = ([topology.houses.cul_de_sac_id] == csac);
    topology_csac = topology.houses(houses_csac);
    num_houses = sum(houses_csac);
    house_ids = sort([topology.houses(houses_csac).id]);
    
    % Ground truth
    true_offset = zeros([size(topology_csac), 1]);
    ground_truth = table([topology_csac.id]', true_offset, ...
        [topology_csac.service_pipe_insulation_W_m_K]', ...
        'VariableNames', {'house_id', 'true_offset', 'true_U'});
    
    % Meter data for this CSAC
    csac_table = table([topology_csac.id]', [topology_csac.service_pipe_len_m]', ...
        [topology_csac.dist_on_cul_de_sac_m]', ...
        'VariableNames', {'house_id', 'length_service_m', 'x_pos_m'});
    meter_data_csac = meter_data(ismember(meter_data.house_id, house_ids), :);
    meter_data_csac = join(meter_data_csac, csac_table, "Keys", "house_id");
    
    % Initialize filters
    ukf_states = cell(1, num_houses);
    pf_particles = cell(1, num_houses);
    pf_states = cell(1, num_houses);
    ukf_states_snapshot = cell(1, num_houses);
    pf_states_snapshot = cell(1, num_houses);
    
    for i = 1:num_houses
        x_init = [randn() * 0.3; 0.12 + randn() * 0.03];
        ukf_states{i} = initialize_kalman_filter(x_init, params.Q_base, params.R_base, params.P_base);
        pf_particles{i} = initialize_pf_state(params.num_particles, x_init, params.P_base, params.config);
        pf_states{i}.x = x_init;
        pf_states{i}.P = params.P_base;
        ukf_states_snapshot{i}.x = x_init;
        ukf_states_snapshot{i}.P = params.P_base;
        pf_states_snapshot{i}.x = x_init;
        pf_states_snapshot{i}.P = params.P_base;
    end
    
    % Pack everything into csac_state
    csac_state.U_csac = U_csac;
    csac_state.num_houses = num_houses;
    csac_state.house_ids = house_ids;
    csac_state.ground_truth = ground_truth;
    csac_state.meter_data_csac = meter_data_csac;
    
    csac_state.ukf_states = ukf_states;
    csac_state.pf_particles = pf_particles;
    csac_state.pf_states = pf_states;
    csac_state.ukf_states_snapshot = ukf_states_snapshot;
    csac_state.pf_states_snapshot = pf_states_snapshot;
    csac_state.snapshot_timestep = NaN;
    
    csac_state.ukf_offsets = nan(num_houses, 1);
    csac_state.pf_offsets = nan(num_houses, 1);
    
    csac_state.last_valid_T_main_ukf_C = nan(num_houses, 1);
    csac_state.last_valid_T_main_pf_C = nan(num_houses, 1);
    csac_state.last_update_timestamp_ukf = NaT(num_houses, 1);
    csac_state.last_update_timestamp_pf = NaT(num_houses, 1);
    
    % Gate statistics
    csac_state.gate_reject_count_ukf = 0;
    csac_state.gate_reject_count_pf = 0;
    csac_state.gate_accept_count_ukf = 0;
    csac_state.gate_accept_count_pf = 0;
    
    % Season tracking
    csac_state.season_active = false;
    csac_state.was_season_active = false;
    csac_state.season_end_timestep = NaN;
    csac_state.season_start_timestep = NaN;
    csac_state.season_start_date = NaT;
    csac_state.season_count = 0;
    csac_state.last_season_end_date = NaT;
    csac_state.last_season_duration_days = 0;
    csac_state.season_gate_active_days = 0;
    csac_state.season_gate_inactive_days = 0;
    csac_state.last_gate_check_date = NaT;
    csac_state.cumulative_winter_days = 0;
    csac_state.last_dormant_P_growth_date = NaT;
    
    % Loggers
    csac_state.logger_ukf = initialize_logger(num_houses, length(timestamps), house_ids);
    csac_state.logger_pf = initialize_logger(num_houses, length(timestamps), house_ids);
end