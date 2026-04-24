function [pass, fail, details] = test_generate_true_trajectories()
    pass = 0; fail = 0; details = {};

    gt = table([1;2], [0.5; -0.3], [0.12; 0.10], ...
        'VariableNames', {'house_id', 'true_offset', 'true_U'});
    timestamps = datetime(2019,1,1) + hours(0:8759); % 1 year hourly

    %% Test 1: 'none' gives constant trajectories
    try
        dc.type = 'none';
        dc.offset_drift_per_year = 0;
        dc.step_time = NaT;
        dc.offset_step = 0;
        traj = generate_true_trajectories(gt, timestamps, dc);
        assert_true(numel(traj) == 2, 'should have 2 houses');
        assert_near(traj{1}.offset(1), 0.5, 1e-10, 'house 1 offset start');
        assert_near(traj{1}.offset(end), 0.5, 1e-10, 'house 1 offset end');
        assert_near(traj{2}.U(1), 0.10, 1e-10, 'house 2 U constant');
        assert_near(traj{2}.U(end), 0.10, 1e-10, 'house 2 U end');
        pass = pass + 1;
        fprintf('  ✓ None gives constant trajectories\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: 'linear' drift increases offset over time
    try
        dc.type = 'linear';
        dc.offset_drift_per_year = 0.5;
        dc.step_time = NaT;
        dc.offset_step = 0;
        traj = generate_true_trajectories(gt, timestamps, dc);
        % After 1 year, offset should increase by ~0.5
        assert_near(traj{1}.offset(end), 0.5 + 0.5, 0.01, 'house 1 offset after 1 year');
        assert_near(traj{1}.offset(1), 0.5, 1e-10, 'house 1 offset at start');
        % U should be unchanged
        assert_near(traj{1}.U(end), 0.12, 1e-10, 'U should not drift');
        pass = pass + 1;
        fprintf('  ✓ Linear drift increases offset\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: 'step' changes offset at step_time
    try
        dc.type = 'step';
        dc.offset_drift_per_year = 0;
        dc.step_time = datetime(2019,7,1);
        dc.offset_step = 0.8;
        traj = generate_true_trajectories(gt, timestamps, dc);
        % Before step
        idx_before = find(timestamps < dc.step_time, 1, 'last');
        assert_near(traj{1}.offset(idx_before), 0.5, 1e-10, 'before step should be base');
        % After step
        idx_after = find(timestamps >= dc.step_time, 1, 'first');
        assert_near(traj{1}.offset(idx_after), 0.5 + 0.8, 1e-10, 'after step should be base+step');
        pass = pass + 1;
        fprintf('  ✓ Step changes offset at correct time\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: U-value never drifts regardless of type
    try
        dc.type = 'linear';
        dc.offset_drift_per_year = 1.0;
        dc.step_time = NaT;
        dc.offset_step = 0;
        traj = generate_true_trajectories(gt, timestamps, dc);
        assert_near(traj{1}.U(1), traj{1}.U(end), 1e-10, 'U should never change');
        pass = pass + 1;
        fprintf('  ✓ U-value constant regardless of drift type\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end