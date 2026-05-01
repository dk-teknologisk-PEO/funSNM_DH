function cs = make_csac_with_offsets(offsets, positions, P_offset_val)
%MAKE_CSAC_WITH_OFFSETS Creates a minimal CSAC state struct for testing.
%
%   Args:
%       offsets (1xN double): Offset values per house.
%       positions (1xN double): Position on CSAC per house [m].
%       P_offset_val (scalar): Covariance value for offset (small = converged).

    num_houses = numel(offsets);
    cs.num_houses = num_houses;
    cs.house_ids = 0:(num_houses - 1);

    cs.ukf_states = cell(1, num_houses);
    for i = 1:num_houses
        cs.ukf_states{i}.x = [offsets(i); 0.12];
        cs.ukf_states{i}.P = diag([P_offset_val, 0.001]);
    end

    % Create minimal meter_data with positions
    house_id_col = cs.house_ids(:);
    x_pos_col = positions(:);
    ts = repmat(datetime(2019,1,1), num_houses, 1);
    flow = repmat(200, num_houses, 1);
    T_supply = repmat(68, num_houses, 1);

    cs.meter_data = table(house_id_col, ts, flow, T_supply, x_pos_col, ...
        'VariableNames', {'house_id', 'timestamp', 'flow_kg_h', 'T_supply_C', 'x_pos_m'});

    cs.gate_accept_count = 100;
    cs.T_inlet_fitted = NaN;
    cs.current_total_flow = sum(flow);
end