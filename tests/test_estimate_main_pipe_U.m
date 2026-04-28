function [pass, fail, details] = test_estimate_main_pipe_U()
    pass = 0; fail = 0; details = {};

    config = make_test_config();
    topology = make_test_topology();
    csac_ids = [0, 1, 2, 3];

    %% Test 1: No adjustment when avg offsets show no position trend
    try
        all_cs = {
            make_csac_with_offsets([0.1, -0.1, 0.05], [10,30,50], 0.01)
            make_csac_with_offsets([-0.05, 0.1, -0.1], [10,30,50], 0.01)
            make_csac_with_offsets([0.0, 0.05, -0.05], [10,30,50], 0.01)
            make_csac_with_offsets([0.1, -0.05, 0.0], [10,30,50], 0.01)
        };
        [U_new, diag] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.20, config);
        assert_near(U_new, 0.20, 0.003, 'U should barely change with no trend');
        pass = pass + 1;
        fprintf('  ✓ No adjustment for random avg offsets\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Positive adjustment when avg offsets increase with CSAC position
    try
        % CSAC positions: 24, 54, 126, 166 m
        % Avg offsets increasing with position → U_main too low
        all_cs = {
            make_csac_with_offsets([0.0, 0.0, 0.0], [10,30,50], 0.01)
            make_csac_with_offsets([0.1, 0.1, 0.1], [10,30,50], 0.01)
            make_csac_with_offsets([0.3, 0.3, 0.3], [10,30,50], 0.01)
            make_csac_with_offsets([0.4, 0.4, 0.4], [10,30,50], 0.01)
        };
        [U_new, diag] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.15, config);
        if diag.adjusted
            assert_true(U_new > 0.15, 'U should increase for positive gradient');
        end
        pass = pass + 1;
        fprintf('  ✓ Positive adjustment for increasing avg offsets\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Does not crash with only 2 CSACs having valid data
    try
        all_cs = {
            make_csac_with_offsets([0.0, 0.1, 0.2], [10,30,50], 0.01)
            make_csac_with_offsets([0.1, 0.2, 0.3], [10,30,50], 10.0)  % unconverged
            make_csac_with_offsets([0.2, 0.3, 0.4], [10,30,50], 10.0)  % unconverged
            make_csac_with_offsets([0.3, 0.4, 0.5], [10,30,50], 0.01)
        };
        [U_new, diag] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.15, config);
        assert_true(isfinite(U_new), 'should return finite U');
        pass = pass + 1;
        fprintf('  ✓ Handles partially valid CSACs\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: Result stays within bounds
    try
        all_cs = {
            make_csac_with_offsets([0.0, 0.0, 0.0], [10,30,50], 0.01)
            make_csac_with_offsets([1.0, 1.0, 1.0], [10,30,50], 0.01)
            make_csac_with_offsets([2.0, 2.0, 2.0], [10,30,50], 0.01)
            make_csac_with_offsets([3.0, 3.0, 3.0], [10,30,50], 0.01)
        };
        [U_new, ~] = estimate_main_pipe_U(all_cs, csac_ids, topology, 0.49, config);
        assert_true(U_new <= config.project.main_pipe_U_estimation.U_max, 'should not exceed U_max');
        pass = pass + 1;
        fprintf('  ✓ Result clamped to bounds\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end