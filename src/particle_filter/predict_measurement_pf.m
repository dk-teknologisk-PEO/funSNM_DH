% predict_measurement_pf.m

function [z_pred, P_zz, propagated_particles] = predict_measurement_pf(particles, house_data, T_ambient_C, R, Q, config)
%PREDICT_MEASUREMENT_PF Predicts measurement and its variance using the particle cloud.
%
%   Performs the prediction/propagation step of the PF and then calculates the
%   statistics of the resulting measurement distribution.

    N = size(particles, 2); % Number of particles

    U_min = config.project.bounds.U_min;
    U_max = config.project.bounds.U_max;
    offset_min = config.project.bounds.offset_min;
    offset_max = config.project.bounds.offset_max;


    % 1. Prediction (Propagate each particle)
    % Each particle evolves according to the process model (a random walk).
    process_noise = chol(Q)' * randn(2, N);
    propagated_particles = particles + process_noise;

    % Apply State Constraints to the propagated particles
    % offset_min = -2.0; offset_max = 2.0;
    % U_min = 0.10; U_max = 0.20;
    propagated_particles(1, :) = max(min(propagated_particles(1, :), offset_max), offset_min);
    propagated_particles(2, :) = max(min(propagated_particles(2, :), U_max), U_min);

    % 2. Predict Measurement for Each Particle
    predicted_measurements = zeros(1, N);
    for i = 1:N
        particle_offset = propagated_particles(1, i);
        particle_U      = propagated_particles(2, i);
        
        predicted_measurements(i) = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h, ...
            particle_U, house_data.length_service_m, T_ambient_C) - particle_offset;
    end

    % 3. Calculate Prediction Statistics
    z_pred = mean(predicted_measurements);
    P_zz_from_particles = var(predicted_measurements);

    % 4. Calculate Total Innovation Variance
    P_zz = P_zz_from_particles + R;
end