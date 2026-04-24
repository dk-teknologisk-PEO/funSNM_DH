function [pass, fail, details] = test_apply_offset_drift_to_data()
    pass = 0; fail = 0; details = {};

    timestamps = datetime(2019,1,1) + hours(0:8759);

    %% Test 1: 'none' returns unchanged data
    try
        data = make_simple_meter_table(timestamps);
        T_before = data.T_supply_C;
        dc.type = 'none';
        dc.offset_drift_per_year = 0;
        dc.step_time = NaT;
        dc.offset_step = 0;
        data = apply_offset_drift_to_data(data, [1], dc, timestamps);
        assert_near(max(abs(data.T_supply_C - T_before)), 0, 1e-10, 'data should not change');
        pass = pass + 1;
        fprintf('  ✓ None returns unchanged data\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: 'linear' drift reduces T_supply over time
    try
        data = make_simple_meter_table(timestamps);
        T_before_start = data.T_supply_C(1);
        T_before_end = data.T_supply_C(end);
        dc.type = 'linear';
        dc.offset_drift_per_year = 0.5;
        dc.step_time = NaT;
        dc.offset_step = 0;
        data = apply_offset_drift_to_data(data, [1], dc, timestamps);
        % Start should be unchanged (t_years = 0)
        assert_near(data.T_supply_C(1), T_before_start, 1e-6, 'start should be unchanged');
        % End should be reduced by ~0.5 (1 year of drift)
        expected_reduction = 0.5 * years(timestamps(end) - timestamps(1));
        actual_reduction = T_before_end - data.T_supply_C(end);
        assert_near(actual_reduction, expected_reduction, 0.01, 'end should be reduced by drift amount');
        pass = pass + 1;
        fprintf('  ✓ Linear drift reduces T_supply\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: 'step' changes T_supply after step_time
    try
        data = make_simple_meter_table(timestamps);
        step_time = datetime(2019,7,1);
        dc.type = 'step';
        dc.offset_drift_per_year = 0;
        dc.step_time = step_time;
        dc.offset_step = 1.0;
        T_before = data.T_supply_C;
        data = apply_offset_drift_to_data(data, [1], dc, timestamps);

        before_mask = data.timestamp < step_time;
        after_mask = data.timestamp >= step_time;
        % Before step: unchanged
        assert_near(max(abs(data.T_supply_C(before_mask) - T_before(before_mask))), 0, 1e-10, ...
            'before step should be unchanged');
        % After step: reduced by 1.0
        diffs = T_before(after_mask) - data.T_supply_C(after_mask);
        assert_near(mean(diffs), 1.0, 1e-10, 'after step should be reduced by 1.0');
        pass = pass + 1;
        fprintf('  ✓ Step changes T_supply after step_time\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: Linear drift magnitude is correct at midpoint
    try
        data = make_simple_meter_table(timestamps);
        dc.type = 'linear';
        dc.offset_drift_per_year = 1.0; % 1°C per year
        dc.step_time = NaT;
        dc.offset_step = 0;
        T_before = data.T_supply_C;
        data = apply_offset_drift_to_data(data, [1], dc, timestamps);

        mid_idx = round(numel(timestamps) / 2);
        t_years_mid = years(timestamps(mid_idx) - timestamps(1));
        expected_reduction = 1.0 * t_years_mid;
        actual_reduction = T_before(mid_idx) - data.T_supply_C(mid_idx);
        assert_near(actual_reduction, expected_reduction, 0.01, 'midpoint reduction should match');
        pass = pass + 1;
        fprintf('  ✓ Linear drift magnitude correct at midpoint\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end