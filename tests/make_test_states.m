function [states, snapshots] = make_test_states(n)
%MAKE_TEST_STATES Creates n test UKF state structs and snapshots.
    states = cell(1, n);
    snapshots = cell(1, n);
    for i = 1:n
        states{i}.x = [0.1*i; 0.12];
        states{i}.P = diag([0.5, 0.01]);
        states{i}.Q = diag([1e-6, 1e-8]);
        states{i}.R = 0.01;
        snapshots{i}.x = states{i}.x;
        snapshots{i}.P = states{i}.P;
    end
end