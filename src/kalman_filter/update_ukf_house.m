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
    
    % UKF parameters
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
    
    sensitivity_offset = 1.0;
    
    dU = 0.001;
    temp_nominal = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, current_U, house_data.length_service_m, T_soil_C);
    temp_perturbed = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, current_U + dU, house_data.length_service_m, T_soil_C);
    sensitivity_U = abs((temp_perturbed - temp_nominal) / dU);
    
    std_offset = sqrt(P(1,1));
    std_U      = sqrt(P(2,2));
    
    impact_offset = sensitivity_offset * std_offset;
    impact_U      = sensitivity_U * std_U;
    
    epsilon = 1e-9;
    focus_on_U = impact_U / (impact_U + impact_offset + epsilon);
    flow_quality = max(0, min(1, (alpha_flow - alpha_min) / (1 - alpha_min)));
    
    focus_on_U = focus_on_U * flow_quality;
    focus_on_offset = 1 - focus_on_U;

    noise_scaler = 10;
    
    Q_offset_var = Q_base_offset_var * (1 + noise_scaler * focus_on_offset);
    Q_U_var      = Q_base_U_var      * (1 + noise_scaler * focus_on_U);
    
    Q = diag([Q_offset_var, Q_U_var]);    

    %% 3. UKF Prediction Step
    % Process model is identity — no mean reversion
    x_pred = x;

    P_pred = P + Q;

    alpha_forget = config.project.initialization.ukf.alpha_forget;
    P_max_offset = 4.0^2;
    P_max_U      = 0.2^2;

    P_pred(1,1) = min(alpha_forget * P_pred(1,1), P_max_offset);
    P_pred(2,2) = min(alpha_forget * P_pred(2,2), P_max_U);

    % Covariance floor
    P_floor_offset = (0.05)^2;
    P_floor_U      = (0.005)^2;

    P_pred(1,1) = max(P_pred(1,1), P_floor_offset);
    P_pred(2,2) = max(P_pred(2,2), P_floor_U);

    % Ensure P_pred is symmetric and positive definite before sigma points
    P_pred = (P_pred + P_pred') / 2;
    P_pred = ensure_positive_definite(P_pred);

    %% 4. UKF Measurement Update Step
    [sigma_points, Wm, Wc] = generate_sigma_points(x_pred, P_pred, alpha, beta, kappa);
    
    num_sigma = size(sigma_points, 2);
    Z_sigma = zeros(1, num_sigma);
    for i = 1:num_sigma
        sigma_offset = sigma_points(1, i);
        sigma_U      = sigma_points(2, i);
        
        T_service_out = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, sigma_U, house_data.length_service_m, T_soil_C);
        Z_sigma(i) = T_service_out - sigma_offset;
    end
    
    z_pred = sum(Wm .* Z_sigma, 2);
    
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
    
    %% 5. State and Covariance Update
    K = P_xz / P_zz;
    
    y = house_data.T_supply_C - z_pred;
    
    x_new = x_pred + K * y;

    % Maximum step size limiter
    max_dx = [0.15; 0.01];
    dx = x_new - x_pred;
    was_limited = false;
    for idx = 1:length(dx)
        if abs(dx(idx)) > max_dx(idx)
            dx(idx) = sign(dx(idx)) * max_dx(idx);
            was_limited = true;
        end
    end
    x_new = x_pred + dx;

    % Joseph form covariance update (numerically stable)
    % P_new = (I - K*H) * P_pred * (I - K*H)' + K * R * K'
    % For UKF, the effective H is implicit. We use the equivalent form:
    % P_new = P_pred - K * P_zz * K'  (standard form)
    % then enforce positive definiteness.
    P_new = P_pred - K * P_zz * K';
    P_new = (P_new + P_new') / 2;
    P_new = ensure_positive_definite(P_new);

    % If step was limited, inflate covariance
    if was_limited
        P_new(1,1) = P_new(1,1) * 1.5;
        P_new(2,2) = P_new(2,2) * 1.5;
    end
    
    % Apply physical constraints
    offset_min = config.project.cutoff.offset_min;
    offset_max = config.project.cutoff.offset_max;
    U_min = config.project.cutoff.U_min;
    U_max = config.project.cutoff.U_max;
    
    x_new(1) = max(min(x_new(1), offset_max), offset_min);
    x_new(2) = max(min(x_new(2), U_max), U_min);

    %% 6. Pack Results
    state.x = x_new;
    state.P = P_new;
    state.y = y;

    diagnostics.y = y;
    diagnostics.P_zz = P_zz;
    diagnostics.K = K;
    diagnostics.P_pred_diag = diag(P_pred);
    diagnostics.P_post_diag = diag(P_new);
end


function P = ensure_positive_definite(P)
%ENSURE_POSITIVE_DEFINITE Forces a 2x2 matrix to be positive definite.
%   Uses eigenvalue decomposition. Cheap for 2x2.
    min_eig = 1e-10;  % minimum eigenvalue
    [V, D] = eig(P);
    d = diag(D);
    if any(d < min_eig)
        d = max(d, min_eig);
        P = V * diag(d) * V';
        P = (P + P') / 2;  % ensure symmetry after reconstruction
    end
end