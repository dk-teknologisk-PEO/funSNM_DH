function [pass, fail, details] = test_apply_master_offset()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    %% Test 1: No offset applied when within deadzone
    try
        [states, ~] = make_test_states(3);
        x_before = states{1}.x(1);
        states = apply_master_offset(states, 0.001, 5, config);
        assert_near(states{1}.x(1), x_before, 1e-10, 'offset should not change in deadzone');
        pass = pass + 1;
        fprintf('  ✓ Deadzone prevents small offset\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Offset applied when above deadzone
    try
        [states, ~] = make_test_states(3);
        x_before = states{1}.x(1);
        states = apply_master_offset(states, 1.0, 5, config);
        assert_true(states{1}.x(1) ~= x_before, 'offset should change above deadzone');
        pass = pass + 1;
        fprintf('  ✓ Offset applied above deadzone\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Not applied when too few active houses
    try
        [states, ~] = make_test_states(3);
        x_before = states{1}.x(1);
        states = apply_master_offset(states, 1.0, 1, config); % only 1 active
        assert_near(states{1}.x(1), x_before, 1e-10, 'should not apply with too few houses');
        pass = pass + 1;
        fprintf('  ✓ Not applied with insufficient active houses\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: NaN offset does nothing
    try
        [states, ~] = make_test_states(3);
        x_before = states{1}.x(1);
        states = apply_master_offset(states, NaN, 5, config);
        assert_near(states{1}.x(1), x_before, 1e-10, 'NaN offset should do nothing');
        pass = pass + 1;
        fprintf('  ✓ NaN offset handled\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: Offset is clamped to max
    try
        [states, ~] = make_test_states(3);
        x_before = states{1}.x(1);
        gamma = config.project.master_offset.gamma;
        % Apply a huge offset — should be clamped to 0.5
        states = apply_master_offset(states, 100.0, 5, config);
        expected_shift = gamma * 0.5; % clamped to 0.5 then damped
        assert_near(states{1}.x(1), x_before - expected_shift, 1e-10, 'offset should be clamped');
        pass = pass + 1;
        fprintf('  ✓ Large offset clamped correctly\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 6: All houses get the same shift
    try
        [states, ~] = make_test_states(5);
        x_before = zeros(5,1);
        for i = 1:5, x_before(i) = states{i}.x(1); end
        states = apply_master_offset(states, 0.3, 5, config);
        shifts = zeros(5,1);
        for i = 1:5, shifts(i) = states{i}.x(1) - x_before(i); end
        for i = 2:5
            assert_near(shifts(i), shifts(1), 1e-10, 'all houses should get same shift');
        end
        pass = pass + 1;
        fprintf('  ✓ All houses get identical shift\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end