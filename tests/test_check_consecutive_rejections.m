function [pass, fail, details] = test_check_consecutive_rejections()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    %% Test 1: Counter increments on rejection
    try
        state = make_test_ukf_state();
        counter = 0;
        [~, counter] = check_consecutive_rejections(state, counter, true, false, config);
        assert_true(counter == 1, 'counter should be 1 after one rejection');
        [~, counter] = check_consecutive_rejections(state, counter, true, false, config);
        assert_true(counter == 2, 'counter should be 2 after two rejections');
        pass = pass + 1;
        fprintf('  ✓ Counter increments on rejection\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Counter resets on acceptance
    try
        state = make_test_ukf_state();
        counter = 5;
        [~, counter] = check_consecutive_rejections(state, counter, false, true, config);
        assert_true(counter == 0, 'counter should reset on acceptance');
        pass = pass + 1;
        fprintf('  ✓ Counter resets on acceptance\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: P inflated after max consecutive rejections
    try
        state = make_test_ukf_state();
        state.P = diag([0.01, 0.001]); % small P
        P_before = state.P;
        counter = 9; % one away from threshold (max_consecutive=10)
        [state, counter] = check_consecutive_rejections(state, counter, true, false, config);
        assert_true(state.P(1,1) > P_before(1,1), 'P(1,1) should be inflated');
        assert_true(state.P(2,2) > P_before(2,2), 'P(2,2) should be inflated');
        assert_true(counter == 0, 'counter should reset after inflation');
        pass = pass + 1;
        fprintf('  ✓ P inflated after max consecutive rejections\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: P capped at maximum
    try
        state = make_test_ukf_state();
        state.P = diag([3.0, 0.05]); % already large
        counter = 9;
        [state, counter] = check_consecutive_rejections(state, counter, true, false, config);
        P_max_offset = config.project.consecutive_rejection.P_max_offset^2;
        P_max_U = config.project.consecutive_rejection.P_max_U^2;
        assert_true(state.P(1,1) <= P_max_offset + 1e-10, 'P(1,1) should be capped');
        assert_true(state.P(2,2) <= P_max_U + 1e-10, 'P(2,2) should be capped');
        pass = pass + 1;
        fprintf('  ✓ P capped at maximum\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: Skip (neither accepted nor rejected) leaves counter unchanged
    try
        state = make_test_ukf_state();
        counter = 3;
        [~, counter] = check_consecutive_rejections(state, counter, false, false, config);
        assert_true(counter == 3, 'counter should not change on skip');
        pass = pass + 1;
        fprintf('  ✓ Skip leaves counter unchanged\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 6: No inflation before reaching threshold
    try
        state = make_test_ukf_state();
        state.P = diag([0.01, 0.001]);
        P_before = state.P;
        counter = 0;
        for k = 1:9 % one less than threshold
            [state, counter] = check_consecutive_rejections(state, counter, true, false, config);
        end
        assert_near(state.P(1,1), P_before(1,1), 1e-10, 'P should not change before threshold');
        assert_true(counter == 9, 'counter should be 9');
        pass = pass + 1;
        fprintf('  ✓ No inflation before threshold\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end