function gp = make_test_gate_params()
%MAKE_TEST_GATE_PARAMS Creates default gate params for testing.
    gp.last_valid_T_main_C = 70.0;
    gp.max_delta_T_change = 5.0;
    gp.hard_innovation_gate_C = 2.0;
    gp.max_nis = 5.0;
    gp.disable_stability_gate = false;
    gp.debug_print = false;
    gp.house_id = 1;
    gp.time = datetime(2019,11,15,12,0,0);
end