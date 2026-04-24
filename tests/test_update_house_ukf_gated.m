function [pass, fail, details] = test_update_house_ukf_gated()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    %% Test 1: Stability gate blocks update when T_main changes too fast
    try
        state = make_test_ukf_state();
        house_data = make_test_house_data();
        gate_params = make_test_gate_params();
        gate_params.last_valid_T_main_C = 50.0; % Big jump from house_data.T_main_ukf_C (~70)
        gate_params.max_delta_T_change = 5.0;
        gate_params.disable_stability_gate = false;

        [~, result] = update_house_ukf_gated(state, house_data, 8.0, config, gate_params);
        assert_true(result.skipped, 'should be skipped by stability gate');
        assert_true(strcmp(result.reason, 'stability_gate'), 'reason should be stability_gate');
        assert_false(result.accepted, 'should not be accepted');
        assert_false(result.rejected, 'should not be rejected');
        pass = pass + 1;
        fprintf('  ✓ Stability gate blocks large T_main change\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Stability gate skipped when disabled
    try
        state = make_test_ukf_state();
        house_data = make_test_house_data();
        gate_params = make_test_gate_params();
        gate_params.last_valid_T_main_C = 50.0;
        gate_params.max_delta_T_change = 5.0;
        gate_params.disable_stability_gate = true;

        [~, result] = update_house_ukf_gated(state, house_data, 8.0, config, gate_params);
        assert_false(result.skipped, 'should not be skipped when gate disabled');
        pass = pass + 1;
        fprintf('  ✓ Stability gate disabled correctly\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Hard innovation gate rejects large innovation
    try
        state = make_test_ukf_state();
        house_data = make_test_house_data();
        house_data.T_supply_C = 50.0; % Way off from predicted
        gate_params = make_test_gate_params();
        gate_params.hard_innovation_gate_C = 2.0;

        [~, result] = update_house_ukf_gated(state, house_data, 8.0, config, gate_params);
        assert_true(result.rejected, 'should be rejected by hard innovation gate');
        assert_true(strcmp(result.reason, 'hard_innovation_gate'), 'reason should be hard_innovation_gate');
        pass = pass + 1;
        fprintf('  ✓ Hard innovation gate rejects large innovation\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: Normal update is accepted
    try
        state = make_test_ukf_state();
        house_data = make_test_house_data();
        gate_params = make_test_gate_params();

        x_before = state.x;
        [state_new, result] = update_house_ukf_gated(state, house_data, 8.0, config, gate_params);
        assert_true(result.accepted, 'normal update should be accepted');
        assert_true(strcmp(result.reason, 'accepted'), 'reason should be accepted');
        assert_true(~isempty(fieldnames(result.diagnostics)), 'diagnostics should be populated');
        pass = pass + 1;
        fprintf('  ✓ Normal update accepted\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: NaN stability gate — first update with no previous T_main
    try
        state = make_test_ukf_state();
        house_data = make_test_house_data();
        gate_params = make_test_gate_params();
        gate_params.last_valid_T_main_C = NaN;

        [~, result] = update_house_ukf_gated(state, house_data, 8.0, config, gate_params);
        assert_false(result.skipped, 'NaN last_valid should not trigger stability gate');
        pass = pass + 1;
        fprintf('  ✓ NaN last_valid_T_main handled correctly\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 6: State unchanged after rejection
    try
        state = make_test_ukf_state();
        house_data = make_test_house_data();
        house_data.T_supply_C = 50.0; % Force rejection
        gate_params = make_test_gate_params();
        gate_params.hard_innovation_gate_C = 2.0;

        x_before = state.x;
        P_before = state.P;
        [state_after, ~] = update_house_ukf_gated(state, house_data, 8.0, config, gate_params);
        assert_near(state_after.x(1), x_before(1), 1e-10, 'x(1) should not change after rejection');
        assert_near(state_after.x(2), x_before(2), 1e-10, 'x(2) should not change after rejection');
        assert_near(state_after.P(1,1), P_before(1,1), 1e-10, 'P should not change after rejection');
        pass = pass + 1;
        fprintf('  ✓ State unchanged after rejection\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end