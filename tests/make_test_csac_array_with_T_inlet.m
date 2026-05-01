function all_cs = make_test_csac_array_with_T_inlet(T_inlets)
%MAKE_TEST_CSAC_ARRAY_WITH_T_INLET Creates a cell array of minimal CSAC
%   structs with specified T_inlet_fitted values.

    num_csacs = numel(T_inlets);
    all_cs = cell(1, num_csacs);

    for c = 1:num_csacs
        cs = make_csac_with_offsets([0, 0, 0], [10, 30, 50], 0.01);
        cs.T_inlet_fitted = T_inlets(c);
        cs.current_total_flow = 500;
        cs.gate_accept_count = 100;
        cs.season_state = initialize_season_state();
        cs.season_state.active = isfinite(T_inlets(c));
        all_cs{c} = cs;
    end
end