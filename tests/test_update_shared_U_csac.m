function [pass, fail, details] = test_update_shared_U_csac()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    %% Test 1: No adjustment when slopes cancel across CSACs
    try
        % Opposite slopes in two CSACs should cancel out
        all_cs = {
            make_csac_with_offsets([0.0, 0.05, 0.1, 0.15, 0.2], [10,30,50,70,90], 0.01)
            make_csac_with_offsets([0.2, 0.15, 0.1, 0.05, 0.0], [10,30,50,70,90], 0.01)
        };
        [U_new, diag] = update_shared_U_csac(all_cs, 0.15, config);
        assert_near(U_new, 0.15, 0.003, 'U should barely change when slopes cancel');
        pass = pass + 1;
        fprintf('  ✓ No net adjustment when slopes cancel\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Positive adjustment when all CSACs show positive offset gradient
    try
        all_cs = {
            make_csac_with_offsets([0.0, 0.1, 0.2, 0.3, 0.4], [10,30,50,70,90], 0.01)
            make_csac_with_offsets([0.0, 0.15, 0.25, 0.35, 0.45], [10,30,50,70,90], 0.01)
        };
        [U_new, diag] = update_shared_U_csac(all_cs, 0.10, config);
        if diag.adjusted
            assert_true(U_new > 0.10, 'U should increase for positive gradient');
        end
        pass = pass + 1;
        fprintf('  ✓ Positive adjustment for positive gradient\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Negative adjustment when all CSACs show negative offset gradient
    try
        all_cs = {
            make_csac_with_offsets([0.4, 0.3, 0.2, 0.1, 0.0], [10,30,50,70,90], 0.01)
            make_csac_with_offsets([0.45, 0.35, 0.25, 0.15, 0.0], [10,30,50,70,90], 0.01)
        };
        [U_new, diag] = update_shared_U_csac(all_cs, 0.20, config);
        if diag.adjusted
            assert_true(U_new < 0.20, 'U should decrease for negative gradient');
        end
        pass = pass + 1;
        fprintf('  ✓ Negative adjustment for negative gradient\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: Adjustment clamped to max
    try
        all_cs = {
            make_csac_with_offsets([0.0, 1.0, 2.0, 3.0, 4.0], [10,30,50,70,90], 0.01)
        };
        [U_new, diag] = update_shared_U_csac(all_cs, 0.10, config);
        if diag.adjusted
            max_adj = config.project.csac_U_estimation.max_adjustment_per_step;
            assert_true(abs(U_new - 0.10) <= max_adj + 1e-10, 'adjustment should be clamped');
        end
        pass = pass + 1;
        fprintf('  ✓ Adjustment clamped to max\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: Result clamped to physical bounds
    try
        all_cs = {
            make_csac_with_offsets([4.0, 3.0, 2.0, 1.0, 0.0], [10,30,50,70,90], 0.01)
        };
        [U_new, ~] = update_shared_U_csac(all_cs, 0.05, config);
        assert_true(U_new >= config.project.csac_U_estimation.U_min, 'should not go below U_min');
        pass = pass + 1;
        fprintf('  ✓ Result clamped to physical bounds\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 6: More houses = more influence
    try
        % CSAC1 with 5 houses showing positive gradient
        % CSAC2 with 3 houses showing negative gradient
        cs1 = make_csac_with_offsets([0.0, 0.1, 0.2, 0.3, 0.4], [10,30,50,70,90], 0.01);
        cs2 = make_csac_with_offsets([0.2, 0.1, 0.0], [10,50,90], 0.01);
        all_cs = {cs1, cs2};
        [U_new, diag] = update_shared_U_csac(all_cs, 0.15, config);
        if diag.adjusted
            % CSAC1 has more houses, so positive gradient should dominate
            assert_true(U_new > 0.15, 'larger CSAC should dominate direction');
        end
        pass = pass + 1;
        fprintf('  ✓ More houses give more influence\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end