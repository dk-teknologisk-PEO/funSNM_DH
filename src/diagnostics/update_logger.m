function logger = update_logger(logger, t_idx, house_idx, time, state, diagnostics)
% UPDATE_LOGGER Adds data from a single update step to the logger.
    
    logger.timestamps(t_idx) = time;
    
    % Store state and posterior covariance
    logger.state_estimates(:, house_idx, t_idx) = state.x;
    logger.covariance_posterior(:, house_idx, t_idx) = diag(state.P);
    
    % Store diagnostics
    logger.innovations(house_idx, t_idx) = diagnostics.y;
    logger.innovation_variance(house_idx, t_idx) = diagnostics.P_zz;
    logger.nis(house_idx, t_idx) = diagnostics.y' * (diagnostics.P_zz \ diagnostics.y);
    logger.kalman_gain(:, house_idx, t_idx) = diagnostics.K;
end