% In initialize_kalman.m
function state = initialize_kalman_filter(x_init, Q_base, R_base, P_base)
    state.x = x_init; % [offset; U_value]
    state.P = P_base;% diag([2^2, 0.2^2]);
    state.Q = Q_base;% diag([0.001^2, 0.0001^2]); % This will be adapted later
    state.R = R_base;% (config.project.initialization.measurement_noise).^2;
    state.y = NaN;
end