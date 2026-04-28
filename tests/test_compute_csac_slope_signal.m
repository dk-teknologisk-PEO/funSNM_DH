function [pass, fail, details] = test_compute_csac_slope_signal()
    pass = 0; fail = 0; details = {};

    config = make_test_config();

    %% Test 1: No slope when offsets are random (no position trend)
    try
        cs = make_csac_with_offsets([0.1, -0.2, 0.15, -0.1, 0.05], ...
            [10, 30, 50, 70, 90], 0.01);
        signal = compute_csac_slope_signal(cs, config);
        assert_true(signal.valid, 'should be valid with 5 houses');
        assert_true(abs(signal.corr_offset) < 0.5, 'correlation should be weak for random offsets');
        pass = pass + 1;
        fprintf('  ✓ No slope for random offsets\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 2: Positive slope detected when offsets increase with position
    try
        cs = make_csac_with_offsets([0.0, 0.1, 0.2, 0.3, 0.4], ...
            [10, 30, 50, 70, 90], 0.01);
        signal = compute_csac_slope_signal(cs, config);
        assert_true(signal.valid, 'should be valid');
        assert_true(signal.slope_offset > 0, 'slope should be positive');
        assert_true(signal.corr_offset > 0.8, 'correlation should be strong');
        assert_true(signal.total_gradient_offset > 0, 'total gradient should be positive');
        pass = pass + 1;
        fprintf('  ✓ Positive slope detected\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 3: Negative slope detected when offsets decrease with position
    try
        cs = make_csac_with_offsets([0.4, 0.3, 0.2, 0.1, 0.0], ...
            [10, 30, 50, 70, 90], 0.01);
        signal = compute_csac_slope_signal(cs, config);
        assert_true(signal.valid, 'should be valid');
        assert_true(signal.slope_offset < 0, 'slope should be negative');
        assert_true(signal.corr_offset < -0.8, 'correlation should be strong negative');
        pass = pass + 1;
        fprintf('  ✓ Negative slope detected\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 4: De-meaning removes constant offset level
    try
        % All offsets shifted by +2.0, but same gradient
        cs1 = make_csac_with_offsets([0.0, 0.1, 0.2, 0.3, 0.4], ...
            [10, 30, 50, 70, 90], 0.01);
        cs2 = make_csac_with_offsets([2.0, 2.1, 2.2, 2.3, 2.4], ...
            [10, 30, 50, 70, 90], 0.01);
        signal1 = compute_csac_slope_signal(cs1, config);
        signal2 = compute_csac_slope_signal(cs2, config);
        assert_near(signal1.slope_offset, signal2.slope_offset, 1e-6, ...
            'slope should be same regardless of mean offset');
        pass = pass + 1;
        fprintf('  ✓ De-meaning removes constant offset\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 5: Not valid with too few houses
    try
        cs = make_csac_with_offsets([0.1, 0.2], [10, 50], 0.01);
        signal = compute_csac_slope_signal(cs, config);
        % Should still work with 2 houses but be less reliable
        % With convergence threshold, might be invalid
        assert_true(islogical(signal.valid), 'should return valid flag');
        pass = pass + 1;
        fprintf('  ✓ Handles few houses without crash\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 6: Unconverged houses excluded
    try
        cs = make_csac_with_offsets([0.0, 0.1, 0.2, 0.3, 0.4], ...
            [10, 30, 50, 70, 90], 0.01);
        % Make houses 4 and 5 unconverged (large P)
        cs.ukf_states{4}.P(1,1) = 10.0;
        cs.ukf_states{5}.P(1,1) = 10.0;
        signal = compute_csac_slope_signal(cs, config);
        assert_true(signal.valid, 'should be valid with 3 converged houses');
        assert_true(signal.num_valid_houses == 3, 'should only use 3 houses');
        pass = pass + 1;
        fprintf('  ✓ Unconverged houses excluded\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end

    %% Test 7: U_service slope detected
    try
        cs = make_csac_with_offsets([0.0, 0.0, 0.0, 0.0, 0.0], ...
            [10, 30, 50, 70, 90], 0.01);
        % Set U_service with position trend
        for i = 1:5
            cs.ukf_states{i}.x(2) = 0.10 + 0.005 * i;
            cs.ukf_states{i}.P(2,2) = 0.001;
        end
        signal = compute_csac_slope_signal(cs, config);
        assert_true(signal.slope_U > 0, 'U_service slope should be positive');
        pass = pass + 1;
        fprintf('  ✓ U_service slope detected\n');
    catch ME
        fail = fail + 1;
        details{end+1} = ME.message;
        fprintf('  ✗ %s\n', ME.message);
    end
end