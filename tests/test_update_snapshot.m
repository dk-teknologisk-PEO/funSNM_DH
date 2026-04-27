function [pass, fail, details] = test_update_snapshot()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    %% Test 1: No snapshot saved before min season duration
    try
        s = initialize_season_state();
        s.start_date = datetime(2019,11,1);
        current_date = datetime(2019,11,10); % only 9 days
        [states, ~] = make_test_states(3);
        snaps_before = cell(size(states));
        for i = 1:3
            snaps_before{i}.x = [99; 99]; % sentinel
            snaps_before{i}.P = diag([99, 99]);
        end
        [snaps, ts] = update_snapshot(s, current_date, 100, states, snaps_before, config);
        assert_near(snaps{1}.x(1), 99, 1e-10, 'snapshot should not be updated');
        assert_true(isnan(ts), 'snapshot_timestep should remain NaN');
        pass = pass + 1;
        fprintf('  ✓ No snapshot before min season duration\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Snapshot saved after min season duration
    try
        s = initialize_season_state();
        s.start_date = datetime(2019,10,1);
        current_date = datetime(2019,12,15); % ~75 days
        [states, ~] = make_test_states(3);
        states{1}.x = [0.5; 0.13];
        snaps_before = cell(size(states));
        for i = 1:3
            snaps_before{i}.x = [0; 0];
            snaps_before{i}.P = diag([1, 1]);
        end
        [snaps, ts] = update_snapshot(s, current_date, 200, states, snaps_before, config);
        assert_near(snaps{1}.x(1), 0.5, 1e-10, 'snapshot should match current state');
        assert_true(ts == 200, 'snapshot_timestep should be updated');
        pass = pass + 1;
        fprintf('  ✓ Snapshot saved after min duration\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end