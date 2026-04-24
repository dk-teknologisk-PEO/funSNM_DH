function state = make_test_ukf_state()
%MAKE_TEST_UKF_STATE Creates a single test UKF state struct.
    state.x = [0.0; 0.12];
    state.P = diag([1.0, 0.01]);
    state.Q = diag([1e-6, 1e-8]);
    state.R = 0.01;
end