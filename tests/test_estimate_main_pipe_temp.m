function [pass, fail, details] = test_estimate_main_pipe_temp()
    pass = 0; fail = 0; details = {};

    config = make_test_config();
    topology = make_test_topology();

    %% Test 1: Returns valid result with enough active CSACs
    try
        all_cs = make_test_csac_array_with_T_inlet([72, 71, 70, 69]);
        csac_ids = [0, 1, 2, 3];
        [T_inlet, diag] = estimate_main_pipe_temp(all_cs, csac_ids, topology, 8.0, 0.20, config);
        assert_true(diag.valid, 'should be valid with 4 CSACs');
        assert_true(all(isfinite(T_inlet)), 'all T_inlet should be finite');
        pass = pass + 1;
        fprintf('  ✓ Valid result with 4 CSACs\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Temperatures decrease with position
    try
        all_cs = make_test_csac_array_with_T_inlet([72, 71, 70, 69]);
        csac_ids = [0, 1, 2, 3];
        [T_inlet, ~] = estimate_main_pipe_temp(all_cs, csac_ids, topology, 8.0, 0.20, config);
        for c = 2:4
            assert_true(T_inlet(c) <= T_inlet(c-1) + 0.1, ...
                'temperature should generally decrease with position');
        end
        pass = pass + 1;
        fprintf('  ✓ Temperatures decrease with position\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Network inlet temperature is higher than all CSAC inlets
    try
        all_cs = make_test_csac_array_with_T_inlet([72, 71, 70, 69]);
        csac_ids = [0, 1, 2, 3];
        [T_inlet, diag] = estimate_main_pipe_temp(all_cs, csac_ids, topology, 8.0, 0.20, config);
        assert_true(diag.T_network_inlet >= max(T_inlet), ...
            'network inlet should be >= all CSAC inlets');
        pass = pass + 1;
        fprintf('  ✓ Network inlet higher than CSAC inlets\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: Invalid with too few active CSACs
    try
        all_cs = make_test_csac_array_with_T_inlet([72, NaN, NaN, NaN]);
        csac_ids = [0, 1, 2, 3];
        [~, diag] = estimate_main_pipe_temp(all_cs, csac_ids, topology, 8.0, 0.20, config);
        assert_false(diag.valid, 'should be invalid with only 1 active CSAC');
        pass = pass + 1;
        fprintf('  ✓ Invalid with too few CSACs\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: Higher U_main gives more temperature drop
    try
        all_cs = make_test_csac_array_with_T_inlet([72, 71, 70, 69]);
        csac_ids = [0, 1, 2, 3];
        [T_low_U, ~] = estimate_main_pipe_temp(all_cs, csac_ids, topology, 8.0, 0.10, config);
        [T_high_U, ~] = estimate_main_pipe_temp(all_cs, csac_ids, topology, 8.0, 0.30, config);
        % With higher U, the fitted network inlet should be higher to compensate
        % But the temperature drop between first and last CSAC should be larger
        drop_low = T_low_U(1) - T_low_U(4);
        drop_high = T_high_U(1) - T_high_U(4);
        assert_true(drop_high > drop_low, 'higher U should give more temperature drop');
        pass = pass + 1;
        fprintf('  ✓ Higher U_main gives more temperature drop\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end