function [pass, fail, details] = test_initialize_season_state()
    pass = 0; fail = 0; details = {};

    %% Test 1: All fields exist and have correct defaults
    try
        s = initialize_season_state();
        assert_false(s.active, 'active should be false');
        assert_true(s.count == 0, 'count should be 0');
        assert_true(isnat(s.start_date), 'start_date should be NaT');
        assert_true(isnan(s.start_timestep), 'start_timestep should be NaN');
        assert_true(isnan(s.end_timestep), 'end_timestep should be NaN');
        assert_true(isnat(s.last_gate_check_date), 'last_gate_check_date should be NaT');
        assert_true(isnat(s.last_active_date), 'last_active_date should be NaT');
        assert_true(isnat(s.last_real_season_end_date), 'last_real_season_end_date should be NaT');
        assert_true(s.last_season_duration_days == 0, 'last_season_duration_days should be 0');
        assert_false(s.cooldown_active, 'cooldown_active should be false');
        assert_true(s.active_days == 0, 'active_days should be 0');
        assert_true(s.inactive_days == 0, 'inactive_days should be 0');
        assert_true(isnan(s.snapshot_timestep), 'snapshot_timestep should be NaN');
        pass = pass + 1;
        fprintf('  ✓ All fields initialized correctly\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Two calls produce independent structs
    try
        s1 = initialize_season_state();
        s2 = initialize_season_state();
        s1.active = true;
        s1.count = 5;
        assert_false(s2.active, 'modifying s1 should not affect s2');
        assert_true(s2.count == 0, 's2.count should still be 0');
        pass = pass + 1;
        fprintf('  ✓ Independent struct instances\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end