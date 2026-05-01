% run_all_tests.m
% Runs all unit tests and reports results.

clear all; close all; clc;

addpath('../src/kalman_filter', '../src/network_model', '../src/data_handling', ...
    '../src/diagnostics', '../src/gates', '../src/CSACs', '../config')

fprintf('========================================\n');
fprintf('RUNNING UNIT TESTS\n');
fprintf('========================================\n\n');

total_pass = 0;
total_fail = 0;
total_skip = 0;
failed_tests = {};

test_functions = {
    @test_initialize_season_state
    @test_manage_heating_season
    @test_apply_season_actions
    @test_update_house_ukf_gated
    @test_apply_master_offset
    @test_update_snapshot
    @test_compute_house_kpis
    @test_generate_true_trajectories
    @test_apply_offset_drift_to_data
    @test_build_daily_T_air_max_table
    @test_print_csac_summary
    @test_check_consecutive_rejections
    @test_compute_csac_slope_signal
    @test_update_shared_U_csac
    @test_estimate_main_pipe_temp
    @test_estimate_main_pipe_U
};

for i = 1:numel(test_functions)
    func = test_functions{i};
    func_name = func2str(func);
    fprintf('--- %s ---\n', func_name);
    try
        [pass, fail, details] = func();
        total_pass = total_pass + pass;
        total_fail = total_fail + fail;
        if fail > 0
            for d = 1:numel(details)
                failed_tests{end+1} = sprintf('%s: %s', func_name, details{d}); %#ok<AGROW>
            end
        end
    catch ME
        total_fail = total_fail + 1;
        failed_tests{end+1} = sprintf('%s: CRASHED — %s', func_name, ME.message); %#ok<AGROW>
        fprintf('  CRASHED: %s\n', ME.message);
    end
    fprintf('\n');
end

fprintf('========================================\n');
fprintf('RESULTS: %d passed, %d failed\n', total_pass, total_fail);
fprintf('========================================\n');

if ~isempty(failed_tests)
    fprintf('\nFAILED TESTS:\n');
    for i = 1:numel(failed_tests)
        fprintf('  ✗ %s\n', failed_tests{i});
    end
end

if total_fail == 0
    fprintf('\nAll tests passed.\n');
end