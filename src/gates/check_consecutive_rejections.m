function [ukf_state, rejection_counter] = check_consecutive_rejections(ukf_state, rejection_counter, was_rejected, was_accepted, config)
%CHECK_CONSECUTIVE_REJECTIONS Inflates covariance after repeated gate rejections.
%
%   If a house is repeatedly rejected by the innovation or NIS gate while
%   preconditions (flow, stability, etc.) are met, it signals that the
%   state estimate is wrong rather than the measurements being bad.
%   After enough consecutive rejections, the covariance is inflated to
%   allow the filter to re-converge.
%
%   Args:
%       ukf_state (struct): Current UKF state for this house.
%       rejection_counter (int): Number of consecutive rejections so far.
%       was_rejected (logical): Whether the current update was rejected.
%       was_accepted (logical): Whether the current update was accepted.
%       config (struct): Project configuration.
%
%   Returns:
%       ukf_state (struct): Possibly modified UKF state (inflated P).
%       rejection_counter (int): Updated counter.

    reject_cfg = config.project.consecutive_rejection;
    max_consecutive = reject_cfg.max_consecutive;
    P_inflation_factor = reject_cfg.P_inflation_factor;
    P_max_offset = reject_cfg.P_max_offset^2;  % stored as std, square for variance
    P_max_U = reject_cfg.P_max_U^2;

    if was_accepted
        % Reset counter on any successful update
        rejection_counter = 0;
        return;
    end

    if was_rejected
        rejection_counter = rejection_counter + 1;

        if rejection_counter >= max_consecutive
            % Inflate covariance to allow re-convergence
            ukf_state.P(1,1) = min(ukf_state.P(1,1) * P_inflation_factor, P_max_offset);
            ukf_state.P(2,2) = min(ukf_state.P(2,2) * P_inflation_factor, P_max_U);

            fprintf('  CONSEC REJECT INFLATE | %d consecutive rejections | P_offset=%.4f, P_U=%.6f\n', ...
                rejection_counter, sqrt(ukf_state.P(1,1)), sqrt(ukf_state.P(2,2)));

            % Reset counter so it can build up again if still rejected
            rejection_counter = 0;
        end
    end
    % If neither rejected nor accepted (skipped), counter unchanged
end