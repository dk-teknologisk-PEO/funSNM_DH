function [particles, state_estimate, state_covariance, diagnostics] = update_pf_house(particles, house_data, T_ambient_C, R, Q)
%UPDATE_PF_HOUSE Performs a particle filter update for a single house.
%   Uses a Sequential Importance Resampling (SIR) algorithm.
%
%   Args:
%       particles (2xN matrix): The current set of N particles [offset; U_value].
%       house_data (table row): Current meter data.
%       T_main_C (scalar): Estimated main pipe temperature.
%       T_ambient_C (scalar): Ambient (soil) temperature.
%       R (scalar): Measurement noise variance.
%       Q (2x2 matrix): Process noise covariance.
%
%   Returns:
%       particles (2xN matrix): The updated set of particles.
%       state_estimate (2x1 vector): The new mean state estimate.
%       state_covariance (2x2 matrix): The new state covariance.

    N = size(particles, 2); % Number of particles
    weights = ones(1, N) / N; % Initialize weights for this step

    %% 1. Prediction (Propagate each particle)
    % Each particle evolves according to the process model (a random walk).
    % We add noise sampled from the process noise covariance Q.
    process_noise = chol(Q)' * randn(2, N);
    particles = particles + process_noise;

    %% 2. Apply State Constraints (Project particles back into valid range)
    % This is crucial for preventing particles from wandering into absurd regions.
    offset_min = -2.0; offset_max = 2.0;
    U_min = 0.10; U_max = 0.20;
    
    particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
    particles(2, :) = max(min(particles(2, :), U_max), U_min);
    
    %% 3. Weighting (Evaluate each particle against the measurement)
    % Calculate the likelihood of the measurement given each particle's state.
    
    measurement = house_data.T_supply_C;
    likelihoods = zeros(1, N);
    
    for i = 1:N
        % Get the state hypothesis from the current particle
        particle_offset = particles(1, i);
        particle_U      = particles(2, i);
        
        % Predict the measurement using the particle's state (REUSING YOUR FUNCTION)
        predicted_measurement = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, particle_U, house_data.length_service_m, T_ambient_C) - particle_offset;
        
        % Calculate the probability of the actual measurement under a Gaussian assumption
        % This is the core of the weighting step.
        error = measurement - predicted_measurement;
        likelihoods(i) = (1 / sqrt(2 * pi * R)) * exp(-0.5 * error^2 / R);
    end
    
    % Update the weights by multiplying with the new likelihoods
    weights = weights .* likelihoods;
    
    % Normalize the weights so they sum to 1
    sum_weights = sum(weights);
    if sum_weights > 1e-15 % Avoid division by zero
        weights = weights / sum_weights;
    else
        % All particles are very unlikely. Re-initialize to avoid collapse.
        weights = ones(1, N) / N;
    end

    %% 4. Resampling (Based on effective sample size)
    % This step combats particle degeneracy.
    
    N_eff = 1 / sum(weights.^2);
    resample_threshold = N * 0.5; % Resample if N_eff drops below 50%
    
    if N_eff < resample_threshold
        % Perform systematic resampling (a low-variance method)
        indices = zeros(1, N);
        C = cumsum(weights);
        u1 = rand() / N;
        i = 1;
        for j = 1:N
            u = u1 + (j-1)/N;
            while u > C(i)
                i = i + 1;
            end
            indices(j) = i;
        end
        
        % Replace old particles with resampled ones and reset weights
        particles = particles(:, indices);
        weights = ones(1, N) / N;
        
        % OPTIONAL BUT RECOMMENDED: Jittering / Regularization
        % Add a small amount of noise to the resampled particles to re-introduce diversity.
        jitter_Q = Q * 0.1; % Use a fraction of the process noise
        particles = particles + chol(jitter_Q)' * randn(2, N);
        % Re-apply constraints after jittering
        particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
        particles(2, :) = max(min(particles(2, :), U_max), U_min);
    end

    %% 5. State Estimation (Compute mean and covariance from particles)
    % The final estimate is the weighted mean of the particles.
    state_estimate = particles * weights';
    
    % The covariance is the weighted covariance of the particle cloud.
    state_covariance = zeros(2, 2);
    for i = 1:N
        diff = particles(:, i) - state_estimate;
        state_covariance = state_covariance + weights(i) * (diff * diff');
    end
    %% 6. Calculate Diagnostics for Logging
    
    % a) Predict the measurement based on the FINAL state estimate
    predicted_measurement = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, state_estimate(2), house_data.length_service_m, T_ambient_C) - state_estimate(1);
    
    % b) Calculate the residual (innovation) 'y'
    diagnostics.y = house_data.T_supply_C - predicted_measurement;
    
    % c) Estimate the innovation variance 'P_zz'
    % The uncertainty in the measurement prediction comes from two sources:
    % 1. The uncertainty in the state estimate (propagated through the measurement function).
    % 2. The measurement noise R.
    % For simplicity, we can approximate this. A simple but effective way is to
    % just use the measurement noise itself, as the particle cloud already
    % accounts for state uncertainty. A more complex method would linearize
    % around the mean, but that defeats the purpose of a PF.
    % So, we'll use a pragmatic approximation.
    diagnostics.P_zz = state_covariance(1,1) + R; % Simplified: uncertainty from offset + measurement noise
    
    % d) Populate other fields with NaNs since they don't exist in a PF
    diagnostics.K = [NaN; NaN];
    diagnostics.P_pred_diag = [NaN; NaN];
    diagnostics.P_post_diag = diag(state_covariance); % We do have the posterior P
end

