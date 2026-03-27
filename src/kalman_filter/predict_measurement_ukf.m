% predict_measurement_ukf.m

function [z_pred, P_zz] = predict_measurement_ukf(state, house_data, T_soil_C)
%PREDICT_MEASUREMENT_UKF Predicts the measurement and its variance without updating the state.
%   This function performs the first half of the UKF update step.

    % Unpack state and UKF parameters
    x_pred = state.x; % Process model is x_k+1 = x_k
    P_pred = state.P + state.Q; % A priori covariance
    alpha = 1e-3; beta = 2; kappa = 0; n = length(x_pred);

    % 1. Generate sigma points from the a priori distribution
    [sigma_points, Wm, Wc] = generate_sigma_points(x_pred, P_pred, alpha, beta, kappa);
    
    % 2. Propagate sigma points through the measurement function h(x)
    num_sigma = size(sigma_points, 2);
    Z_sigma = zeros(1, num_sigma);
    for i = 1:num_sigma
        sigma_offset = sigma_points(1, i);
        sigma_U      = sigma_points(2, i);
        T_service_out = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h, sigma_U, house_data.length_service_m, T_soil_C);
        Z_sigma(i) = T_service_out - sigma_offset;
    end
    
    % 3. Calculate the predicted measurement mean (z_pred)
    z_pred = sum(Wm .* Z_sigma, 2);
    
    % 4. Calculate the innovation covariance (P_zz)
    P_zz_from_state = 0;
    for i = 1:num_sigma
        residual_z = Z_sigma(i) - z_pred;
        P_zz_from_state = P_zz_from_state + Wc(i) * (residual_z * residual_z');
    end
    
    % Add measurement noise to get the final innovation variance
    P_zz = P_zz_from_state + state.R;
end