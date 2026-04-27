function [pass, fail, details] = test_apply_offset_drift_to_data()
    pass = 0; fail = 0; details = {};

    timestamps = datetime(2019,1,1) + hours(0:8759);

    %% Test 1: 'none' returns unchanged data
    try
        data = make_simple_meter_table(timestamps);
        T_before = data.T_supply_C;
        dc.type = 'none';
        dc.house_index = 1;
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

    %% Test 2: 'linear' drift reduces T_supply for target house only
    try
        data = make_two_house_meter_table(timestamps);
        T_before = data.T_supply_C;
        dc.type = 'linear';
        dc.house_index = 1; % only house 1 drifts
        dc.offset_drift_per_year = 0.5;
        dc.step_time = NaT;
        dc.offset_step = 0;
        data = apply_offset_drift_to_data(data, [1; 2], dc, timestamps);

        % House 1 should be changed
        h1_mask = data.house_id == 1;
        h1_end_idx = find(h1_mask, 1, 'last');
        t_years_end = years(data.timestamp(h1_end_idx) - timestamps(1));
        expected_reduction = 0.5 * t_years_end;
        actual_reduction = T_before(h1_end_idx) - data.T_supply_C(h1_end_idx);
        assert_near(actual_reduction, expected_reduction, 0.01, 'house 1 end should be reduced');

        % House 2 should be unchanged
        h2_mask = data.house_id == 2;
        assert_near(max(abs(data.T_supply_C(h2_mask) - T_before(h2_mask))), 0, 1e-10, ...
            'house 2 should be unchanged');
        pass = pass + 1;
        fprintf('  ✓ Linear drift applies only to target house\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: 'step' changes only target house after step_time
    try
        data = make_two_house_meter_table(timestamps);
        step_time = datetime(2019,7,1);
        dc.type = 'step';
        dc.house_index = 2; % only house 2 steps
        dc.offset_drift_per_year = 0;
        dc.step_time = step_time;
        dc.offset_step = 1.0;
        T_before = data.T_supply_C;
        data = apply_offset_drift_to_data(data, [1; 2], dc, timestamps);

        % House 1 unchanged
        h1_mask = data.house_id == 1;
        assert_near(max(abs(data.T_supply_C(h1_mask) - T_before(h1_mask))), 0, 1e-10, ...
            'house 1 should be unchanged');

        % House 2 before step: unchanged
        h2_before = data.house_id == 2 & data.timestamp < step_time;
        assert_near(max(abs(data.T_supply_C(h2_before) - T_before(h2_before))), 0, 1e-10, ...
            'house 2 before step should be unchanged');

        % House 2 after step: reduced by 1.0
        h2_after = data.house_id == 2 & data.timestamp >= step_time;
        diffs = T_before(h2_after) - data.T_supply_C(h2_after);
        assert_near(mean(diffs), 1.0, 1e-10, 'house 2 after step should be reduced by 1.0');
        pass = pass + 1;
        fprintf('  ✓ Step applies only to target house\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end