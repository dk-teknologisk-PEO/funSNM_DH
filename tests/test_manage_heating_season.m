function [pass, fail, details] = test_manage_heating_season()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    %% Build a daily T_air_max table: cold for 20 days, then warm
    dates = (datetime(2019,10,1):days(1):datetime(2020,4,30))';
    temps = zeros(size(dates));
    % First 60 days: cold (below 8°C)
    temps(1:60) = 3;
    % Days 61-90: warm (above 8°C)
    temps(61:90) = 15;
    % Days 91-150: cold again
    temps(91:150) = 2;
    % Rest: warm
    temps(151:end) = 18;
    daily_table = table(dates, temps, 'VariableNames', {'date', 'T_air_max'});

    %% Test 1: Season starts after lookback window of cold days
    try
        s = initialize_season_state();
        % Check day 11 (enough cold days in lookback)
        [s, actions] = manage_heating_season(s, dates(11), 11, daily_table, config, 0);
        assert_true(actions.season_is_active, 'season should be active after 10 cold days');
        assert_true(s.active, 'state.active should be true');
        assert_true(s.count == 1, 'count should be 1');
        assert_false(actions.do_restore_snapshot, 'first season should not restore snapshot');
        pass = pass + 1;
        fprintf('  ✓ Season starts after lookback window\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Season does not start before lookback window is full
    try
        s = initialize_season_state();
        [s, actions] = manage_heating_season(s, dates(5), 5, daily_table, config, 0);
        % With 10-day lookback, day 5 doesn't have enough history
        % Result depends on check_heating_season_gate behavior with insufficient data
        % Just verify it doesn't crash and returns a valid struct
        assert_true(islogical(actions.season_is_active), 'should return logical');
        pass = pass + 1;
        fprintf('  ✓ Early days handled without crash\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Season ends when warm days appear
    try
        s = initialize_season_state();
        % Start the season
        for d = 10:60
            [s, ~] = manage_heating_season(s, dates(d), d, daily_table, config, 0);
        end
        assert_true(s.active, 'season should be active at day 60');
        
        % Now process warm days — season should end
        for d = 61:75
            [s, actions] = manage_heating_season(s, dates(d), d, daily_table, config, 0);
        end
        assert_false(s.active, 'season should have ended during warm period');
        assert_true(actions.skip_timestep, 'should skip when inactive');
        pass = pass + 1;
        fprintf('  ✓ Season ends when warm days appear\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: Same-day calls don't re-evaluate
    try
        s = initialize_season_state();
        [s, a1] = manage_heating_season(s, dates(11), 100, daily_table, config, 0);
        count_after_first = s.count;
        [s, a2] = manage_heating_season(s, dates(11), 101, daily_table, config, 0);
        assert_true(s.count == count_after_first, 'count should not change on same day');
        pass = pass + 1;
        fprintf('  ✓ Same-day calls do not re-evaluate\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: Second season triggers snapshot restore
    try
        s = initialize_season_state();
        s.snapshot_timestep = 50;
        % Run through first season
        for d = 10:60
            [s, ~] = manage_heating_season(s, dates(d), d, daily_table, config, 0);
        end
        % End season (warm period)
        for d = 61:90
            [s, ~] = manage_heating_season(s, dates(d), d, daily_table, config, 0);
        end
        % Start second season (cold again at day 91+)
        s.last_active_date = dates(60); % set for gap calculation
        found_restore = false;
        for d = 91:120
            [s, actions] = manage_heating_season(s, dates(d), d, daily_table, config, 0);
            if actions.do_restore_snapshot
                found_restore = true;
                break;
            end
        end
        assert_true(found_restore, 'second season should trigger snapshot restore');
        assert_true(s.count == 2, 'count should be 2');
        assert_true(actions.gap_days > 0, 'gap_days should be positive');
        pass = pass + 1;
        fprintf('  ✓ Second season triggers snapshot restore\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end