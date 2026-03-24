% In a new file: diagnostic_logger.m

function logger = initialize_logger(num_houses, num_timestamps, house_ids)
% INITIALIZE_LOGGER Creates a structure to store all diagnostic variables.
    
    num_states = 2; % [offset, U-value]
    
    % Pre-allocate arrays with NaNs
    logger.timestamps = NaT(1, num_timestamps);
    logger.house_ids = house_ids;
    
    logger.state_estimates = nan(num_states, num_houses, num_timestamps);
    logger.innovations = nan(num_houses, num_timestamps);
    logger.innovation_variance = nan(num_houses, num_timestamps);
    logger.nis = nan(num_houses, num_timestamps);
    logger.kalman_gain = nan(num_states, num_houses, num_timestamps);
    logger.covariance_posterior = nan(num_states, num_houses, num_timestamps);
    
    disp('Diagnostic logger initialized.');
end

