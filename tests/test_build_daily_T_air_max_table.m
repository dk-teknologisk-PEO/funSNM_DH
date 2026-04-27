function [pass, fail, details] = test_build_daily_T_air_max_table()
    pass = 0; fail = 0; details = {};

    %% Test 1: Correct daily max computed
    try
        hours_vec = (0:47)'; % 2 days of hourly data
        T_air_C = table(datetime(2019,10,1) + hours(hours_vec), ...
            [repmat(5,24,1); repmat(10,24,1)], ...
            'VariableNames', {'time', 'values'});
        result = build_daily_T_air_max_table(T_air_C);
        assert_true(height(result) == 2, 'should have 2 days');
        assert_near(result.T_air_max(1), 5, 1e-10, 'day 1 max should be 5');
        assert_near(result.T_air_max(2), 10, 1e-10, 'day 2 max should be 10');
        pass = pass + 1;
        fprintf('  ✓ Daily max computed correctly\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Intra-day variation picks max
    try
        hours_vec = (0:23)';
        temps = sin(2*pi*hours_vec/24) * 5 + 10; % peaks at ~15
        T_air_C = table(datetime(2019,10,1) + hours(hours_vec), temps, ...
            'VariableNames', {'time', 'values'});
        result = build_daily_T_air_max_table(T_air_C);
        assert_near(result.T_air_max(1), max(temps), 1e-10, 'should pick intra-day max');
        pass = pass + 1;
        fprintf('  ✓ Intra-day max picked correctly\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Output is sorted by date
    try
        % Create data out of chronological order
        t1 = datetime(2019,10,2) + hours(0:23)';
        t2 = datetime(2019,10,1) + hours(0:23)';
        T_air_C = table([t1; t2], [repmat(8,24,1); repmat(3,24,1)], ...
            'VariableNames', {'time', 'values'});
        result = build_daily_T_air_max_table(T_air_C);
        assert_true(result.date(1) < result.date(2), 'output should be sorted');
        pass = pass + 1;
        fprintf('  ✓ Output sorted by date\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end