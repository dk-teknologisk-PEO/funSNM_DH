% In update_pf_house.m

function [particles, state_estimate, state_covariance, diagnostics] = update_pf_house(propagated_particles, house_data, T_ambient_C, R, z_pred, config)
%UPDATE_PF_HOUSE Performs a particle filter update (weighting & resampling).
%   Assumes particles have already been propagated through the process model.

    N = size(propagated_particles, 2); % Number of particles
    weights = ones(1, N) / N; % Initialize weights for this step
    particles = propagated_particles; % Use the propagated particles

    % 1. Weighting (Evaluate each particle against the measurement)
    measurement = house_data.T_supply_C;
    likelihoods = zeros(1, N);
    
    for i = 1:N
        particle_offset = particles(1, i);
        particle_U      = particles(2, i);
        predicted_measurement = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, particle_U, house_data.length_service_m, T_ambient_C) - particle_offset;
        error = measurement - predicted_measurement;
        likelihoods(i) = (1 / sqrt(2 * pi * R)) * exp(-0.5 * error^2 / R);
    end
    
    % 2. Normalize weights
    sum_weights = sum(weights .* likelihoods);
    if sum_weights > 1e-15
        weights = (weights .* likelihoods) / sum_weights;
    else
        weights = ones(1, N) / N; % Avoid collapse, re-initialize weights
    end

    % 3. Resampling (if needed)
    N_eff = 1 / sum(weights.^2);
    if N_eff < N * 0.5
        % ... (Systematic resampling logic) ...
        indices = systematic_resample(weights);
        particles = particles(:, indices);
        weights = ones(1, N) / N;
        
        % Optional Jittering
        jitter_Q = [0.001^2, 0; 0, 0.0001^2] * 0.1; % Example small jitter
        particles = particles + chol(jitter_Q)' * randn(2, N);
        particles(1, :) = max(min(particles(1, :), 2.0), -2.0);
        particles(2, :) = max(min(particles(2, :), 0.20), 0.10);
    end

    % 4. State Estimation
    state_estimate = particles * weights';
    state_covariance = zeros(2, 2);
    for i = 1:N
        diff = particles(:, i) - state_estimate;
        state_covariance = state_covariance + weights(i) * (diff * diff');
    end
    
    % 5. Calculate Diagnostics for Logging
    % The innovation 'y' is the difference between the measurement and the
    % a priori prediction 'z_pred' passed into the function. This is the correct
    % value for statistical tests like NIS.
    diagnostics.y = measurement - z_pred;
    
    % We need to estimate the innovation variance 'P_zz'.
    % This was already calculated outside, but for logging purposes we might
    % want to re-calculate it or, more simply, pass it in. For now, let's
    % use a simple approximation based on the new covariance.
    [~, P_zz_approx, ~] = predict_measurement_pf(particles, house_data, T_ambient_C, R, [0;0], config); % Re-predict to get variance from posterior particles
    diagnostics.P_zz = P_zz_approx; 
    
    diagnostics.K = [NaN; NaN];
end

% Helper function for resampling, can be inside the file or separate
function indices = systematic_resample(weights)
    N = length(weights);
    indices = zeros(1, N);
    C = cumsum(weights);
    u1 = rand() / N;
    i = 1;
    for j = 1:N
        u = u1 + (j - 1) / N;
        while u > C(i)
            i = i + 1;
        end
        indices(j) = i;
    end
end