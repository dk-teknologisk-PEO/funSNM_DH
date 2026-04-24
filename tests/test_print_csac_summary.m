function [pass, fail, details] = test_print_csac_summary()
    pass = 0; fail = 0; details = {};

    %% Test 1: Does not crash with normal inputs
    try
        s = initialize_season_state();
        s.active_days = 100;
        s.inactive_days = 265;
        s.count = 2;
        s.end_timestep = 5000;
        print_csac_summary(0, 500, 50, s);
        pass = pass + 1;
        fprintf('  ✓ Normal inputs do not crash\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Does not crash with zero counts
    try
        s = initialize_season_state();
        print_csac_summary(0, 0, 0, s);
        pass = pass + 1;
        fprintf('  ✓ Zero counts do not crash\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Handles active cooldown
    try
        s = initialize_season_state();
        s.cooldown_active = true;
        s.last_real_season_end_date = datetime(2020,3,15);
        print_csac_summary(0, 100, 10, s);
        pass = pass + 1;
        fprintf('  ✓ Active cooldown printed without crash\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end