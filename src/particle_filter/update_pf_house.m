function [particles, state_estimate, state_covariance, diagnostics] = update_pf_house(particles, house_data, T_ambient_C, R, Q, config)
%UPDATE_PF_HOUSE Performs a particle filter update for a single house.
%   Uses a Sequential Importance Resampling (SIR) algorithm with
%   sensitivity-based process noise scaling.
%
%   Args:
%       particles (2xN matrix): Current particles [offset; U_value].
%       house_data (table row): Current meter data.
%       T_ambient_C (scalar): Ambient (soil) temperature.
%       R (scalar): Measurement noise variance.
%       Q (2x2 matrix): Base process noise covariance.
%       config (struct): Configuration structure.
%
%   Returns:
%       particles (2xN matrix): Updated particles.
%       state_estimate (2x1 vector): New mean state estimate.
%       state_covariance (2x2 matrix): New state covariance.
%       diagnostics (struct): Diagnostic information for logging.

    N = size(particles, 2);
    alpha_min = config.project.cutoff.alpha_min;
    
    % Physical bounds
    offset_min = config.project.cutoff.offset_min;
    offset_max = config.project.cutoff.offset_max;
    U_min = config.project.cutoff.U_min;
    U_max = config.project.cutoff.U_max;

    %% 1. Sensitivity-Based Dynamic Process Noise (matching UKF logic)
    % Use the current particle mean as the linearization point
    current_offset = mean(particles(1,:));
    current_U      = mean(particles(2,:));
    std_offset     = std(particles(1,:));
    std_U          = std(particles(2,:));
    
    c = 4186;
    flow_kg_s = house_data.flow_kg_h / 3600;
    theta = (current_U * house_data.length_service_m) / (2 * c);
    alpha_flow = (flow_kg_s - theta) / (flow_kg_s + theta);
    
    % Sensitivity to offset is always 1.0
    sensitivity_offset = 1.0;
    
    % Numerical sensitivity to U-value
    dU = 0.001;
    temp_nominal  = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
        current_U, house_data.length_service_m, T_ambient_C);
    temp_perturbed = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
        current_U + dU, house_data.length_service_m, T_ambient_C);
    sensitivity_U = abs((temp_perturbed - temp_nominal) / dU);
    
    % Normalize by state uncertainty (use particle spread)
    epsilon = 1e-9;
    impact_offset = sensitivity_offset * max(std_offset, 1e-6);
    impact_U      = sensitivity_U * max(std_U, 1e-6);
    
    focus_on_U = impact_U / (impact_U + impact_offset + epsilon);
    
    % Suppress U-updates when flow is barely valid
    flow_quality = max(0, min(1, (alpha_flow - alpha_min) / (1 - alpha_min)));
    focus_on_U = focus_on_U * flow_quality;
    focus_on_offset = 1 - focus_on_U;
    
    % Scale process noise (same logic as UKF)
    noise_scaler = 100;
    Q_scaled = diag([
        Q(1,1) * (1 + noise_scaler * focus_on_offset), ...
        Q(2,2) * (1 + noise_scaler * focus_on_U)
    ]);

    %% 2. Prediction (Propagate particles)
    % Mean reversion on offset (same as UKF)
    lambda_offset = 0.005;
    particles(1,:) = (1 - lambda_offset) * particles(1,:);
    
    % Add process noise
    process_noise = chol(Q_scaled)' * randn(2, N);
    particles = particles + process_noise;

    % Clamp to physical bounds
    particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
    particles(2, :) = max(min(particles(2, :), U_max), U_min);
    
    %% 3. Weighting (Evaluate likelihood of measurement given each particle)
    measurement = house_data.T_supply_C;
    log_likelihoods = -Inf(1, N);  % Pre-allocate with -Inf (zero likelihood default)
    
    for i = 1:N
        particle_offset = particles(1, i);
        particle_U      = particles(2, i);
        
        % Predict measurement using this particle's state
        predicted_measurement = get_supply_temp(house_data.T_main_pf_C, ...
            house_data.flow_kg_h, particle_U, house_data.length_service_m, ...
            T_ambient_C) - particle_offset;
        
        % Guard against NaN from get_supply_temp
        if ~isfinite(predicted_measurement)
            log_likelihoods(i) = -Inf;  % This particle gets zero weight
            continue;
        end
        
        % Log-likelihood under Gaussian measurement model
        err = measurement - predicted_measurement;
        log_likelihoods(i) = -0.5 * err^2 / R;
    end
    
    % --- Log-sum-exp trick for numerical stability ---
    % Find the max among finite values only
    finite_mask = isfinite(log_likelihoods);
    if ~any(finite_mask)
        % All particles produced NaN predictions — catastrophic.
        % Return uniform weights and current particles without update.
        warning('PF: All particles produced non-finite predictions. Skipping update.');
        state_estimate = mean(particles, 2);
        state_covariance = cov(particles');
        diagnostics.y = NaN;
        diagnostics.P_zz = NaN;
        diagnostics.K = [NaN; NaN];
        diagnostics.P_pred_diag = [NaN; NaN];
        diagnostics.P_post_diag = diag(state_covariance);
        diagnostics.N_eff = N;
        return;
    end
    
    max_log_lik = max(log_likelihoods(finite_mask));
    log_weights = log_likelihoods - max_log_lik;
    log_weights(~finite_mask) = -Inf;  % Keep non-finite particles at zero weight
    weights = exp(log_weights);
    
    % Normalize
    sum_weights = sum(weights);
    if sum_weights > 1e-300
        weights = weights / sum_weights;
    else
        % Complete particle collapse — reinitialize with uniform weights
        weights = ones(1, N) / N;
        warning('PF: Weight collapse detected. Using uniform weights.');
    end

    %% 4. Resampling (Systematic, low-variance)
    N_eff = 1 / sum(weights.^2);
    resample_threshold = N * 0.5;
    
    if N_eff < resample_threshold
        % Systematic resampling
        indices = systematic_resample(weights, N);
        
        particles = particles(:, indices);
        weights = ones(1, N) / N;
        
        % Jittering to re-introduce diversity (use SCALED noise)
        jitter_Q = Q_scaled * 0.1;
        particles = particles + chol(jitter_Q)' * randn(2, N);
        
        % Re-apply constraints after jittering
        particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
        particles(2, :) = max(min(particles(2, :), U_max), U_min);
    end

    %% 5. State Estimation
    state_estimate = particles * weights';
    
    % Weighted covariance
    state_covariance = zeros(2, 2);
    for i = 1:N
        d = particles(:, i) - state_estimate;
        state_covariance = state_covariance + weights(i) * (d * d');
    end
    
    % Ensure minimum variance (prevent particle collapse to a point)
    min_var_offset = 1e-6;
    min_var_U = 1e-8;
    if state_covariance(1,1) < min_var_offset
        state_covariance(1,1) = min_var_offset;
    end
    if state_covariance(2,2) < min_var_U
        state_covariance(2,2) = min_var_U;
    end

    %% 6. Diagnostics — Compute P_zz from particle spread
    % Predict measurement for EACH particle to get proper measurement spread
    predicted_measurements = zeros(1, N);
    for i = 1:N
        predicted_measurements(i) = get_supply_temp(house_data.T_main_pf_C, ...
            house_data.flow_kg_h, particles(2, i), house_data.length_service_m, ...
            T_ambient_C) - particles(1, i);
    end
    
    % Handle any NaN in predictions
    valid_pred = isfinite(predicted_measurements);
    if sum(valid_pred) < 2
        diagnostics.P_zz = R;  % Fallback
    else
        % Use weights for valid particles only (renormalize)
        w_valid = weights(valid_pred);
        w_valid = w_valid / sum(w_valid);
        pred_valid = predicted_measurements(valid_pred);
        
        z_mean = sum(w_valid .* pred_valid);
        P_zz_state = sum(w_valid .* (pred_valid - z_mean).^2);
        diagnostics.P_zz = P_zz_state + R;
    end
    
    % Innovation based on weighted mean prediction
    z_mean_all = sum(weights .* predicted_measurements, 'omitnan');
    diagnostics.y = measurement - z_mean_all;
    diagnostics.K = [NaN; NaN];
    diagnostics.P_pred_diag = [NaN; NaN];
    diagnostics.P_post_diag = diag(state_covariance);
    diagnostics.N_eff = N_eff;
end


%% ========================================================================
% LOCAL FUNCTION: Systematic Resampling
% =========================================================================
function indices = systematic_resample(weights, N)
%SYSTEMATIC_RESAMPLE Low-variance resampling.
    indices = zeros(1, N);
    C = cumsum(weights);
    C(end) = 1.0;  % Ensure exact sum for numerical safety
    
    u1 = rand() / N;
    i = 1;
    for j = 1:N
        u = u1 + (j - 1) / N;
        while u > C(i) && i < N
            i = i + 1;
        end
        indices(j) = i;
    end
end