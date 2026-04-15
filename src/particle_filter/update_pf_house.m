function [particles, state_estimate, state_covariance, diagnostics] = update_pf_house(particles, house_data, T_ambient_C, R, Q_base, config)
%UPDATE_PF_HOUSE Performs a particle filter update for a single house.
%   Uses a Sequential Importance Resampling (SIR) algorithm with
%   sensitivity-based dynamic process noise.
%
%   Args:
%       particles (2xN matrix): The current set of N particles [offset; U_value].
%       house_data (table row): Current meter data.
%       T_ambient_C (scalar): Ambient (soil) temperature.
%       R (scalar): Measurement noise variance.
%       Q_base (2x2 matrix): Base process noise covariance.
%       config (struct): Full configuration struct.
%
%   Returns:
%       particles (2xN matrix): The updated set of particles.
%       state_estimate (2x1 vector): The new mean state estimate.
%       state_covariance (2x2 matrix): The new state covariance.
%       diagnostics (struct): Diagnostic information for logging.

    N = size(particles, 2);
    alpha_min = config.project.cutoff.alpha_min;

    %% 1. Sensitivity-Based Dynamic Process Noise
    % Mirror the UKF's focus factor logic: when flow is high and U is
    % observable, increase U process noise; when flow is low, focus on offset.
    
    current_offset = mean(particles(1, :));
    current_U = mean(particles(2, :));
    current_std_offset = std(particles(1, :));
    current_std_U = std(particles(2, :));
    
    c = 4186;
    flow_kg_s = house_data.flow_kg_h / 3600;
    theta = (current_U * house_data.length_service_m) / (2 * c);
    alpha_flow = (flow_kg_s - theta) / (flow_kg_s + theta);
    
    % Sensitivity of measurement to offset is always 1
    sensitivity_offset = 1.0;
    
    % Numerically approximate sensitivity to U-value
    dU = 0.001;
    temp_nominal = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
        current_U, house_data.length_service_m, T_ambient_C);
    temp_perturbed = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
        current_U + dU, house_data.length_service_m, T_ambient_C);
    sensitivity_U = abs((temp_perturbed - temp_nominal) / dU);
    
    % Normalize by particle cloud spread (analogous to UKF's sqrt(P))
    impact_offset = sensitivity_offset * max(current_std_offset, 0.05);
    impact_U = sensitivity_U * max(current_std_U, 0.005);
    
    epsilon = 1e-9;
    focus_on_U = impact_U / (impact_U + impact_offset + epsilon);
    
    % Scale by flow quality
    flow_quality = max(0, min(1, (alpha_flow - alpha_min) / (1 - alpha_min)));
    focus_on_U = focus_on_U * flow_quality;
    focus_on_offset = 1 - focus_on_U;
    
    % Scale process noise
    noise_scaler = 10;
    Q_offset_var = Q_base(1,1) * (1 + noise_scaler * focus_on_offset);
    Q_U_var = Q_base(2,2) * (1 + noise_scaler * focus_on_U);
    Q = diag([Q_offset_var, Q_U_var]);

    %% 2. Prediction (Propagate particles with mean reversion + dynamic Q)
    lambda_offset = 0.005;
    particles(1, :) = (1 - lambda_offset) * particles(1, :);
    
    process_noise = chol(Q)' * randn(2, N);
    particles = particles + process_noise;

    %% 3. Apply State Constraints
    offset_min = config.project.cutoff.offset_min;
    offset_max = config.project.cutoff.offset_max;
    U_min = config.project.cutoff.U_min;
    U_max = config.project.cutoff.U_max;
    
    particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
    particles(2, :) = max(min(particles(2, :), U_max), U_min);
    
    %% 4. Weighting (Evaluate each particle against measurement)
    measurement = house_data.T_supply_C;
    log_likelihoods = zeros(1, N);
    
    for i = 1:N
        particle_offset = particles(1, i);
        particle_U = particles(2, i);
        
        predicted_measurement = get_supply_temp(house_data.T_main_pf_C, ...
            house_data.flow_kg_h, particle_U, house_data.length_service_m, ...
            T_ambient_C) - particle_offset;
        
        error = measurement - predicted_measurement;
        log_likelihoods(i) = -0.5 * error^2 / R;
    end
    
    % Log-sum-exp trick for numerical stability
    max_log_lik = max(log_likelihoods);
    log_weights = log_likelihoods - max_log_lik;
    weights = exp(log_weights);
    
    sum_weights = sum(weights);
    if sum_weights > 1e-300
        weights = weights / sum_weights;
    else
        weights = ones(1, N) / N;
    end

    %% 5. Resampling (Based on effective sample size)
    N_eff = 1 / sum(weights.^2);
    resample_threshold = N * 0.5;
    
    if N_eff < resample_threshold
        % Systematic resampling
        indices = zeros(1, N);
        C = cumsum(weights);
        C(end) = 1.0;
        u1 = rand() / N;
        i_idx = 1;
        for j = 1:N
            u = u1 + (j - 1) / N;
            while u > C(i_idx)
                i_idx = i_idx + 1;
            end
            indices(j) = i_idx;
        end
        
        particles = particles(:, indices);
        weights = ones(1, N) / N;
        
        % Jittering with dynamic Q (scaled down)
        jitter_Q = Q * 0.1;
        particles = particles + chol(jitter_Q)' * randn(2, N);
        
        % Re-apply constraints
        particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
        particles(2, :) = max(min(particles(2, :), U_max), U_min);
    end

    %% 6. State Estimation with Step Limiting
    % Compute the pre-update mean (before weighting/resampling changed things)
    % This is the predicted mean from step 2, before the measurement update
    pre_update_mean = mean(particles, 2);  % Current cloud center
    
    % Weighted mean (the measurement-updated estimate)
    weighted_mean = particles * weights';
    
    % Step size limiter — if the cloud moved too far, pull it back
    max_dx = [0.15; 0.01];
    dx = weighted_mean - pre_update_mean;
    was_limited = false;
    for idx = 1:length(dx)
        if abs(dx(idx)) > max_dx(idx)
            dx(idx) = sign(dx(idx)) * max_dx(idx);
            was_limited = true;
        end
    end
    
    if was_limited
        % Shift the ENTIRE particle cloud so its weighted mean matches the limited step
        limited_mean = pre_update_mean + dx;
        shift = limited_mean - weighted_mean;
        particles = particles + shift;
        
        % Re-apply constraints after shifting
        particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
        particles(2, :) = max(min(particles(2, :), U_max), U_min);
        
        % Recompute weighted mean from shifted cloud
        state_estimate = particles * weights';
    else
        state_estimate = weighted_mean;
    end
    
    % Weighted covariance
    state_covariance = zeros(2, 2);
    for i = 1:N
        diff = particles(:, i) - state_estimate;
        state_covariance = state_covariance + weights(i) * (diff * diff');
    end

    %% 7. Diagnostics
    z_pred = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
        state_estimate(2), house_data.length_service_m, T_ambient_C) - state_estimate(1);
    
    % Compute P_zz from particle spread in measurement space
    z_particles = zeros(1, N);
    for i = 1:N
        z_particles(i) = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
            particles(2, i), house_data.length_service_m, T_ambient_C) - particles(1, i);
    end
    z_mean = sum(weights .* z_particles);
    P_zz = sum(weights .* (z_particles - z_mean).^2) + R;
    
    diagnostics.y = measurement - z_mean;
    diagnostics.P_zz = P_zz;
    diagnostics.K = [NaN; NaN];
    diagnostics.P_pred_diag = [NaN; NaN];
    diagnostics.P_post_diag = diag(state_covariance);
    diagnostics.N_eff = N_eff;
end