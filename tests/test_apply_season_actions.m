function [pass, fail, details] = test_apply_season_actions()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    P_base = diag([4.0, 0.04]);

    %% Test 1: No action — states unchanged
    try
        [states, snaps] = make_test_states(3);
        actions.do_restore_snapshot = false;
        actions.do_save_snapshot = false;
        actions.gap_days = 0;
        x_before = states{1}.x;
        [states, snaps] = apply_season_actions(actions, states, snaps, P_base, config);
        assert_near(states{1}.x(1), x_before(1), 1e-10, 'state should not change');
        pass = pass + 1;
        fprintf('  ✓ No action leaves states unchanged\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Restore snapshot copies x from snapshot
    try
        [states, snaps] = make_test_states(3);
        % Set snapshot to different values
        for i = 1:3
            snaps{i}.x = [1.0; 0.15];
            snaps{i}.P = diag([0.01, 0.001]);
        end
        actions.do_restore_snapshot = true;
        actions.do_save_snapshot = false;
        actions.gap_days = 0;
        [states, ~] = apply_season_actions(actions, states, snaps, P_base, config);
        assert_near(states{1}.x(1), 1.0, 1e-10, 'x should be restored from snapshot');
        assert_near(states{1}.x(2), 0.15, 1e-10, 'x(2) should be restored from snapshot');
        pass = pass + 1;
        fprintf('  ✓ Restore copies x from snapshot\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Restore with gap inflates P
    try
        [states, snaps] = make_test_states(3);
        for i = 1:3
            snaps{i}.x = [0.5; 0.12];
            snaps{i}.P = diag([0.01, 0.001]);
        end
        actions.do_restore_snapshot = true;
        actions.do_save_snapshot = false;
        actions.gap_days = 100;
        [states, ~] = apply_season_actions(actions, states, snaps, P_base, config);
        assert_true(states{1}.P(1,1) > 0.01, 'P(1,1) should be inflated');
        assert_true(states{1}.P(2,2) > 0.001, 'P(2,2) should be inflated');
        pass = pass + 1;
        fprintf('  ✓ Gap inflates P\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: P never exceeds P_base
    try
        [states, snaps] = make_test_states(3);
        for i = 1:3
            snaps{i}.x = [0; 0.12];
            snaps{i}.P = diag([1.0, 0.01]); % already large
        end
        actions.do_restore_snapshot = true;
        actions.do_save_snapshot = false;
        actions.gap_days = 10000; % huge gap
        [states, ~] = apply_season_actions(actions, states, snaps, P_base, config);
        assert_true(states{1}.P(1,1) <= P_base(1,1) + 1e-10, 'P(1,1) should not exceed P_base');
        assert_true(states{1}.P(2,2) <= P_base(2,2) + 1e-10, 'P(2,2) should not exceed P_base');
        pass = pass + 1;
        fprintf('  ✓ P capped at P_base\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: Save snapshot copies current states
    try
        [states, snaps] = make_test_states(3);
        states{2}.x = [0.77; 0.18];
        states{2}.P = diag([0.05, 0.003]);
        actions.do_restore_snapshot = false;
        actions.do_save_snapshot = true;
        actions.gap_days = 0;
        [~, snaps] = apply_season_actions(actions, states, snaps, P_base, config);
        assert_near(snaps{2}.x(1), 0.77, 1e-10, 'snapshot x should match state');
        assert_near(snaps{2}.P(1,1), 0.05, 1e-10, 'snapshot P should match state');
        pass = pass + 1;
        fprintf('  ✓ Save snapshot copies current states\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end