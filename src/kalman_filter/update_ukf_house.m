% update_ukf_house.m

function [state, diagnostics] = update_ukf_house(state, house_data, T_soil_C, config)
%UPDATE_UKF_HOUSE Performs a UKF update for a single service pipe/meter.

    %% 1. Unpack State and Define Constants
    x = state.x;
    P = state.P;
    Q_base_offset_var = state.Q(1,1);
    Q_base_U_var = state.Q(2,2);
    alpha_min = config.project.cutoff.alpha_min;

    n = length(x);
    
    % UKF tuning parameters
    alpha = 1e-3;
    beta = 2;
    kappa = 0;

    %% 2. Sensitivity-Based Dynamic Process Noise (Q)
    current_offset = x(1);
    current_U      = x(2);

    c = 4186;
    flow_kg_s = house_data.flow_kg_h / 3600;
    theta = (current_U * house_data.length_service_m) / (2 * c);
    alpha_flow = (flow_kg_s - theta) / (flow_kg_s + theta);
    
    % --- Calculate Measurement Sensitivities ---
    sensitivity_offset = 1.0;
    
    dU = 0.001;
    temp_nominal  = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, ...
        current_U, house_data.length_service_m, T_soil_C);
    temp_perturbed = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, ...
        current_U + dU, house_data.length_service_m, T_soil_C);
    sensitivity_U = abs((temp_perturbed - temp_nominal) / dU);
    
    % --- Normalize sensitivities by state uncertainty ---
    std_offset = sqrt(max(P(1,1), 1e-12));
    std_U      = sqrt(max(P(2,2), 1e-12));
    
    impact_offset = sensitivity_offset * std_offset;
    impact_U      = sensitivity_U * std_U;
    
    % --- Define Focus Factor ---
    epsilon = 1e-9;
    focus_on_U = impact_U / (impact_U + impact_offset + epsilon);

    % Suppress U-updates when flow is barely valid
    flow_quality = max(0, min(1, (alpha_flow - alpha_min) / (1 - alpha_min)));
    focus_on_U = focus_on_U * flow_quality;
    focus_on_offset = 1 - focus_on_U;

    % --- Scale Process Noise ---
    % REDUCED from 100 to 10: noise_scaler=100 caused Kalman gains that
    % were far too large, producing offset jumps of 0.3-0.5 per timestep
    % from innovations of only ~0.2. The scaler only needs to create a 
    % contrast between the two states, not massively inflate total noise.
    noise_scaler = 10;
    
    Q_offset_var = Q_base_offset_var * (1 + noise_scaler * focus_on_offset);
    Q_U_var      = Q_base_U_var      * (1 + noise_scaler * focus_on_U);
    
    Q = diag([Q_offset_var, Q_U_var]);    

    %% 3. UKF Prediction Step
    lambda_offset = 0.005;
    x_pred = x;
    x_pred(1) = (1 - lambda_offset) * x(1);

    P_pred = P + Q;

    alpha_forget = config.project.initialization.ukf.alpha_forget;
    P_pred = alpha_forget * P_pred; 

    % --- Ensure P_pred is well-conditioned ---
    P_pred = (P_pred + P_pred') / 2;
    min_var = 1e-10;
    for idx = 1:n
        if P_pred(idx, idx) < min_var
            P_pred(idx, idx) = min_var;
        end
    end

    %% 4. UKF Measurement Update Step
    
    [sigma_points, Wm, Wc] = generate_sigma_points(x_pred, P_pred, alpha, beta, kappa);
    
    num_sigma = size(sigma_points, 2);
    Z_sigma = zeros(1, num_sigma);
    for i = 1:num_sigma
        sigma_offset = sigma_points(1, i);
        sigma_U      = sigma_points(2, i);
        
        % Clamp sigma-point U-values to prevent nonsensical results
        sigma_U = max(sigma_U, config.project.cutoff.U_min * 0.5);
        sigma_U = min(sigma_U, config.project.cutoff.U_max * 2.0);
        
        T_service_out = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, ...
            sigma_U, house_data.length_service_m, T_soil_C);
        
        Z_sigma(i) = T_service_out - sigma_offset;
    end
    
    % Check for NaN in propagated sigma points
    if any(isnan(Z_sigma))
        warning('NaN detected in propagated sigma points. Skipping update.');
        diagnostics.y = NaN;
        diagnostics.P_zz = NaN;
        diagnostics.K = [NaN; NaN];
        diagnostics.P_pred_diag = diag(P_pred);
        diagnostics.P_post_diag = diag(state.P);
        return;
    end
    
    % Predicted measurement mean
    z_pred = sum(Wm .* Z_sigma, 2);
    
    % Innovation covariance (P_zz) and cross-covariance (P_xz)
    P_zz = 0;
    P_xz = zeros(n, 1);
    for i = 1:num_sigma
        residual_z = Z_sigma(i) - z_pred;
        residual_x = sigma_points(:, i) - x_pred;
        
        P_zz = P_zz + Wc(i) * (residual_z * residual_z');
        P_xz = P_xz + Wc(i) * residual_x * residual_z';
    end
    
    % Add measurement noise
    P_zz = P_zz + state.R;
    
    % Guard against degenerate P_zz
    if P_zz < 1e-12
        P_zz = state.R;
    end

    %% 5. State and Covariance Update
    
    % Kalman Gain
    K = P_xz / P_zz;
    K_gated = K;
    
    % --- Soft Kalman gain gating via sigmoid ---
    gate_sharpness = 10; 
    u_update_gain      = 1 / (1 + exp(-gate_sharpness * (focus_on_U - 0.5)));
    offset_update_gain = 1 - u_update_gain;

    K_gated(1) = K(1) * offset_update_gain;
    K_gated(2) = K(2) * u_update_gain;
    
    % Innovation
    y = house_data.T_supply_C - z_pred;
    
    % --- Compute raw state update ---
    dx = K_gated * y;
    
    % --- MAXIMUM STEP SIZE LIMITER ---
    % Prevents any single update from making an excessively large change.
    % This is the last line of defense against outliers that pass the gates.
    max_step_offset = 0.15; % Max offset change per update [°C]
    max_step_U      = 0.01; % Max U-value change per update [W/m/K]
    
    if abs(dx(1)) > max_step_offset
        dx(1) = sign(dx(1)) * max_step_offset;
    end
    if abs(dx(2)) > max_step_U
        dx(2) = sign(dx(2)) * max_step_U;
    end
    
    % Update state with the (possibly limited) step
    x_new = x_pred + dx;
    
    % Covariance update using the original (unlimited) K_gated
    % We still use K_gated (not the limited step) for covariance, because
    % if we limited the step, we want P to remain larger (reflecting that
    % we didn't fully incorporate the measurement information).
    P_new = P_pred - K_gated * P_zz * K_gated';
    P_new = (P_new + P_new') / 2;
    
    % --- Only clamp DIAGONAL elements (variances must be positive) ---
    for idx = 1:n
        if P_new(idx, idx) < 1e-9
            P_new(idx, idx) = 1e-9;
        end
    end
    
    % --- Verify positive-definiteness ---
    [~, flag] = chol(P_new);
    if flag ~= 0
        [V, D] = eig(P_new);
        D = max(D, 1e-9 * eye(n));
        P_new = V * D * V';
        P_new = (P_new + P_new') / 2;
    end
    
    % --- If the step was limited, inflate covariance to reflect that ---
    % --- we didn't fully use the measurement                        ---
    step_was_limited = (abs(K_gated(1) * y) > max_step_offset) || ...
                       (abs(K_gated(2) * y) > max_step_U);
    if step_was_limited
        % Don't reduce P as much — add back some of the innovation info
        P_new = P_new + 0.5 * Q;  % Partial inflation
    end

    % Store state before clamping
    x_unclamped = x_new;

    % --- Physical constraints ---
    offset_min = config.project.cutoff.offset_min;
    offset_max = config.project.cutoff.offset_max;
    U_min = config.project.cutoff.U_min;
    U_max = config.project.cutoff.U_max;
    
    x_new(1) = max(min(x_new(1), offset_max), offset_min);
    x_new(2) = max(min(x_new(2), U_max), U_min);

    offset_is_clamped = (x_new(1) ~= x_unclamped(1));
    U_is_clamped = (x_new(2) ~= x_unclamped(2));

    % --- Covariance inflation at boundaries ---
    if offset_is_clamped
        inflation_factor = 1.1;
        P_new(1,1) = P_new(1,1) * inflation_factor;
        P_new(1,2) = P_new(1,2) * sqrt(inflation_factor);
        P_new(2,1) = P_new(2,1) * sqrt(inflation_factor);
    end
    
    if U_is_clamped
        inflation_factor = 1.1;
        P_new(2,2) = P_new(2,2) * inflation_factor;
        P_new(1,2) = P_new(1,2) * sqrt(inflation_factor);
        P_new(2,1) = P_new(2,1) * sqrt(inflation_factor);
    end

    %% 6. Pack Results
    state.x = x_new;
    state.P = P_new;
    state.y = y;

    diagnostics.y = y;
    diagnostics.P_zz = P_zz;
    diagnostics.K = K_gated;
    diagnostics.P_pred_diag = diag(P_pred);
    diagnostics.P_post_diag = diag(P_new);
end