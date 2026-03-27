% initialize_pf_state.m (Corrected Version)

function particles = initialize_pf_state(N, x_init, P_base)
%INITIALIZE_PF_STATE Creates the initial particle set for the filter.
%   Uses rejection sampling to ensure the initial cloud respects physical
%   constraints while maintaining the shape of the initial Gaussian guess.
%
%   Args:
%       N (integer): Number of particles to generate.
%       x_init (2x1 vector): Initial mean state guess [offset; U_value].
%       P_base (2x2 matrix): Initial covariance of the guess.
%
%   Returns:
%       particles (2xN matrix): The initial cloud of valid particles.

    % Define physical state constraints
    offset_min = -2.0; offset_max = 2.0;
    U_min = 0.10;      U_max = 0.20;

    particles = zeros(2, N);
    num_generated = 0;
    
    % Use Cholesky decomposition for correlated random number generation
    L = chol(P_base, 'lower');

    % Loop until we have generated N valid particles
    while num_generated < N
        % Generate a candidate particle from the Gaussian distribution
        candidate_particle = x_init + L * randn(2, 1);
        
        % Check if the candidate is within all physical bounds
        is_offset_valid = (candidate_particle(1) >= offset_min) && (candidate_particle(1) <= offset_max);
        is_U_valid = (candidate_particle(2) >= U_min) && (candidate_particle(2) <= U_max);
        
        if is_offset_valid && is_U_valid
            % If valid, accept the particle
            num_generated = num_generated + 1;
            particles(:, num_generated) = candidate_particle;
        end
        % If not valid, the loop continues and a new candidate is generated.
    end
end