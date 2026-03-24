function particles = initialize_pf_state(N, x_init, P_base)
%INITIALIZE_PF_STATE Creates the initial particle set for the filter.
%
%   Args:
%       N (integer): Number of particles to generate.
%       x_init (2x1 vector): Initial mean state guess [offset; U_value].
%       P_init (2x2 matrix): Initial covariance of the guess.
%
%   Returns:
%       particles (2xN matrix): The initial cloud of particles.
    
    % Generate N random samples from a Gaussian distribution with the given mean and covariance.
    P_init = P_base; %diag([(config.project.initialization.process_noise_offset)^2, (config.project.initialization.process_noise_U)^2]);
    particles = x_init + chol(P_init)' * randn(2, N);
    
    % Ensure initial particles are within constraints
    offset_min = -2.0; offset_max = 2.0;
    U_min = 0.10; U_max = 0.20;
    particles(1, :) = max(min(particles(1, :), offset_max), offset_min);
    particles(2, :) = max(min(particles(2, :), U_max), U_min);
end