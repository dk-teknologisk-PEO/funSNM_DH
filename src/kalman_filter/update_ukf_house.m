% update_ukf_house.m

function [state, diagnostics] = update_ukf_house(state, house_data, T_soil_C, config)
%UPDATE_UKF_HOUSE Performs a UKF update for a single service pipe/meter.
%   This function estimates the meter temperature offset and the service
%   pipe U-value, using a flow-dependent process noise model.
%
%   Args:
%       state (struct): A struct containing the current filter state:
%                       .x (2x1 vector): [temp_offset; U_value]
%                       .P (2x2 matrix): State covariance
%                       .R (scalar): Measurement noise variance
%       house_data (table row): A single row of meter data containing:
%                               .flow (scalar): Flow in kg/h
%                               .tempHigh (scalar): Measured supply temp
%                               .subPipe (scalar): Service pipe length
%       T_main (scalar): Estimated main pipe temperature at the junction.
%       T_ambient (scalar): Ambient (soil) temperature.
%
%   Returns:
%       state (struct): The updated filter state struct.

    %% 1. Unpack State and Define Constants
    x = state.x; % State: [offset_temp; U_value]
    P = state.P;
    Q_base_offset_var = state.Q(1,1);
    Q_base_U_var = state.Q(2,2);
    alpha_min = config.project.cutoff.alpha_min;

    n = length(x);
    
    % UKF parameters (can be tuned)
    alpha = 1e-3;
    beta = 2;
    kappa = 0;

    %% 2. Sensitivity-Based Dynamic Process Noise (Q)
    
    % Base process noise for slow drift
    % Q_base_offset_var = (0.001)^2;
    % Q_base_U_var      = (0.0001)^2;
    

    % To calculate sensitivity, we need the current best guess for the states
    current_offset = x(1);
    current_U      = x(2);

    c = 4186;
    flow_kg_s = house_data.flow_kg_h / 3600;
    theta = (current_U * house_data.length_service_m) / (2 * c);
    alpha_flow = (flow_kg_s - theta) / (flow_kg_s + theta);
    
    % --- Calculate Measurement Sensitivities ---
    % Sensitivity to offset is always -1. Its magnitude is 1.
    sensitivity_offset = 1.0;
    
    % Numerically approximate the sensitivity to U-value using a small perturbation
    dU = 0.001; % A small change in U-value
    temp_nominal = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, current_U, house_data.length_service_m, T_soil_C);
    temp_perturbed = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, current_U + dU, house_data.length_service_m, T_soil_C);
    sensitivity_U = abs((temp_perturbed - temp_nominal) / dU);
    
    % --- Normalize sensitivities by state uncertainty ---
    % This gives the expected change in measurement per standard deviation of the state.
    std_offset = sqrt(P(1,1));
    std_U      = sqrt(P(2,2));
    
    impact_offset = sensitivity_offset * std_offset;
    impact_U      = sensitivity_U * std_U;
    
    % --- Define Focus Factor ---
    % A value from 0 to 1.
    % focus_on_U = 1 means the measurement is dominated by U-value effects.
    % focus_on_U = 0 means the measurement is dominated by offset effects.
    % Adding a small epsilon to avoid division by zero.
    epsilon = 1e-9;
    focus_on_U = impact_U / (impact_U + impact_offset + epsilon);
    % normalize flow quality above threshold into [0,1]
    flow_quality = max(0, min(1, (alpha_flow - alpha_min) / (1 - alpha_min)));
    
    % if flow is only barely valid, suppress U-updates
    focus_on_U = focus_on_U * flow_quality;
    focus_on_offset = 1 - focus_on_U;


    % --- Scale Process Noise ---
    % If we should focus on the U-value, increase its process noise to allow
    % the filter to update it more aggressively.
    noise_scaler = 10; % Amplifies the effect. Tune this parameter.
    
    Q_offset_var = Q_base_offset_var * (1 + noise_scaler * focus_on_offset);
    Q_U_var      = Q_base_U_var      * (1 + noise_scaler * focus_on_U);
    
    Q = diag([Q_offset_var, Q_U_var]);    



    %% 3. UKF Prediction Step
    % For this problem, the state transition function is the identity,
    % as the parameters are assumed to be near-constant.
    % x_predicted = x;
    % P_predicted = P + Q;
    
    % We generate sigma points directly from the predicted distribution
    % to keep the formulation clean and standard.
    lambda_offset = 0.005; % gentle mean reversion per timestep
    x_pred = x; % Our process model is x_k+1 = x_k
    x_pred(1) = (1-lambda_offset)*x(1);

    P_pred = P + Q;

    alpha_forget = config.project.initialization.ukf.alpha_forget; % Use a small value since updates are sparse (gated)
    P_max_offset = 4.0^2;   % max uncertainty: 4°C std dev
    P_max_U      = 0.2^2;   % max uncertainty: 0.2 W/m/K std dev

    P_pred(1,1) = min(alpha_forget * P_pred(1,1), P_max_offset);
    P_pred(2,2) = min(alpha_forget * P_pred(2,2), P_max_U);

    %% Covariance floor — prevent overconfidence, ensure ongoing responsiveness
    P_floor_offset = (0.05)^2;   % 0.05°C std dev minimum
    P_floor_U      = (0.005)^2;  % 0.005 W/m/K std dev minimum

    P_pred(1,1) = max(P_pred(1,1), P_floor_offset);
    P_pred(2,2) = max(P_pred(2,2), P_floor_U);

    %% 4. UKF Measurement Update Step
    
    % Generate sigma points from the *a priori* (predicted) distribution
    [sigma_points, Wm, Wc] = generate_sigma_points(x_pred, P_pred, alpha, beta, kappa);
    
    % Propagate sigma points through the (nonlinear) measurement function h(x)
    % h(x) = T_supply_at_meter = T_service_out - offset
    num_sigma = size(sigma_points, 2);
    Z_sigma = zeros(1, num_sigma);
    for i = 1:num_sigma
        sigma_offset = sigma_points(1, i);
        sigma_U      = sigma_points(2, i);
        
        % Predict the temperature drop in the service pipe
        T_service_out = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, sigma_U, house_data.length_service_m, T_soil_C);
        
        % Predicted measurement is the pipe outlet temp minus the offset
        Z_sigma(i) = T_service_out - sigma_offset;
    end
    
    % Calculate the predicted measurement mean
    z_pred = sum(Wm .* Z_sigma, 2);
    
    % Calculate innovation covariance (P_zz) and cross-covariance (P_xz)
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
    
    % Calculate Kalman Gain
    K = P_xz / P_zz; % In MATLAB, this is equivalent to P_xz * inv(P_zz)
    
    % Calculate measurement residual (innovation)
    y = house_data.T_supply_C - z_pred;
    
    % Update state and covariance
    x_new = x_pred + K * y;

    %%% STEP 1.2: Maximum step size limiter
    %%% Prevents a single noisy measurement from corrupting a converged estimate.
    max_dx = [0.15; 0.01];  % max change per update: 0.15°C offset, 0.01 W/m/K
    dx = x_new - x_pred;
    was_limited = false;
    for idx = 1:length(dx)
        if abs(dx(idx)) > max_dx(idx)
            dx(idx) = sign(dx(idx)) * max_dx(idx);
            was_limited = true;
        end
    end
    x_new = x_pred + dx;


    P_new = P_pred - K * P_zz * K';
    P_new = (P_new + P_new') / 2; % ensure symmetry

    % Ensure P stays positive definite
    if any(diag(P_new) < 0)
        P_new = P_pred;  % revert to predicted if update broke P
        warning('P went negative — reverting to P_pred for house update');
    end

    %%% STEP 1.2: If step was limited, inflate covariance to signal
    %%% that the filter was overconfident
    if was_limited
        P_new(1,1) = P_new(1,1) * 1.5;
        P_new(2,2) = P_new(2,2) * 1.5;
    end
    
    % Apply physical constraints to the state
    % Define the boundaries for the states
    offset_min = config.project.cutoff.offset_min;%-2.0; 
    offset_max = config.project.cutoff.offset_max;%2.0;
    U_min = config.project.cutoff.U_min;%0.08; 
    U_max = config.project.cutoff.U_max;%0.20;
    
    % Clamp the estimated states to stay within the physical boundaries
    x_new(1) = max(min(x_new(1), offset_max), offset_min); % Clamp offset
    x_new(2) = max(min(x_new(2), U_max), U_min);       % Clamp U-value

    %% 6. Pack Results into State Struct
    state.x = x_new;
    state.P = P_new;
    state.y = y; % Store residual for diagnostics

    diagnostics.y = y;
    diagnostics.P_zz = P_zz;
    diagnostics.K = K;
    diagnostics.P_pred_diag = diag(P_pred);
    diagnostics.P_post_diag = diag(P_new);
end

