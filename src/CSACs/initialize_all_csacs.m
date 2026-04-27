function [all_cs, all_true_traj] = initialize_all_csacs(topology, meter_data, timestamps, ...
    Q_base, R_base, P_base, innovation_gate_initial, config, network_id, drift_config)
%INITIALIZE_ALL_CSACS Initializes state structs for all CSACs in a network.
%
%   Args:
%       topology (struct): Full network topology.
%       meter_data (table): Full meter data for the network.
%       timestamps (datetime array): All unique timestamps.
%       Q_base, R_base, P_base: UKF initialization parameters.
%       innovation_gate_initial (scalar): Initial innovation gate value.
%       config (struct): Project configuration.
%       network_id (int): Network identifier for seeding.
%       drift_config (struct): Drift configuration.
%
%   Returns:
%       all_cs (cell array): One cs struct per CSAC.
%       all_true_traj (cell array): One true_trajectories cell per CSAC.

    csac_ids = [topology.cul_de_sacs.id];
    num_csacs = numel(csac_ids);
    num_timestamps = length(timestamps);

    all_cs = cell(1, num_csacs);
    all_true_traj = cell(1, num_csacs);

    for c = 1:num_csacs
        csac_id = csac_ids(c);
        houses_csac = ([topology.houses.cul_de_sac_id] == csac_id);
        topology_csac = topology.houses(houses_csac);

        % Initialize CSAC state
        cs = initialize_csac_state(topology, topology_csac, houses_csac, ...
            meter_data, num_timestamps, Q_base, R_base, P_base, ...
            innovation_gate_initial, config, network_id);

        % Apply offset drift
        cs.meter_data = apply_offset_drift_to_data(...
            cs.meter_data, cs.house_ids, drift_config, timestamps);

        % Generate true trajectories for KPI evaluation
        true_traj = generate_true_trajectories(...
            cs.ground_truth, timestamps, drift_config);

        all_cs{c} = cs;
        all_true_traj{c} = true_traj;

        fprintf('  Initialized CSAC %d: %d houses\n', csac_id, cs.num_houses);
    end

    fprintf('All %d CSACs initialized.\n\n', num_csacs);
end