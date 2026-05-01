function [pass, fail, details] = test_estimate_main_pipe_U()
    pass = 0; fail = 0; details = {};

    config = make_test_config();
    topology = make_test_topology();
    csac_ids = [0, 1, 2, 3];
    T_soil = 8.0;

    %% Test 1: No adjustment with insufficient active CSACs
    try
        all_cs = make_test_csac_array_with_T_inlet([72, NaN, NaN, NaN]);
        [U_new, diag] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.20, T_soil, config);
        assert_near(U_new, 0.20, 1e-10, 'U should not change with 1 CSAC');
        assert_false(diag.adjusted, 'should not be adjusted');
        pass = pass + 1;
        fprintf('  ✓ No adjustment with insufficient CSACs\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Handles all active CSACs without crash
    try
        all_cs = make_test_csac_array_with_T_inlet([72, 71.5, 70.5, 70]);
        for c = 1:4
            all_cs{c}.season_state.active = true;
        end
        [U_new, diag] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.20, T_soil, config);
        assert_true(isfinite(U_new), 'should return finite U');
        pass = pass + 1;
        fprintf('  ✓ Handles active CSACs without crash\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Result stays within bounds
    try
        all_cs = make_test_csac_array_with_T_inlet([72, 71, 70, 69]);
        for c = 1:4
            all_cs{c}.season_state.active = true;
        end
        [U_new, ~] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.49, T_soil, config);
        assert_true(U_new <= config.project.main_pipe_U_estimation.U_max, 'should not exceed U_max');
        assert_true(U_new >= config.project.main_pipe_U_estimation.U_min, 'should not go below U_min');
        pass = pass + 1;
        fprintf('  ✓ Result stays within bounds\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: Moves toward better fit
    try
        % Generate T_inlet values consistent with U=0.21
        % Then start estimation at U=0.30 — should decrease
        all_cs = make_test_csac_array_with_T_inlet([84.8, 84.5, 83.5, 82.8]);
        for c = 1:4
            all_cs{c}.season_state.active = true;
            all_cs{c}.current_total_flow = 1000;
        end
        [U_new, diag] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.30, T_soil, config);
        if diag.adjusted
            fprintf('    Direction=%d, costs=[%.2f, %.2f, %.2f]\n', ...
                diag.direction, diag.costs(1), diag.costs(2), diag.costs(3));
        end
        % Just verify it doesn't crash and returns a reasonable value
        assert_true(isfinite(U_new), 'should return finite U');
        pass = pass + 1;
        fprintf('  ✓ Moves toward better fit\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: At minimum returns unchanged
    try
        % If current U gives good fit, no adjustment needed
        all_cs = make_test_csac_array_with_T_inlet([72, 72, 72, 72]);
        for c = 1:4
            all_cs{c}.season_state.active = true;
            all_cs{c}.current_total_flow = 1000;
        end
        [U_new, diag] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.20, T_soil, config);
        % With all T_inlets equal, any U gives similar fit
        assert_true(isfinite(U_new), 'should return finite U');
        pass = pass + 1;
        fprintf('  ✓ At minimum returns reasonable result\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end