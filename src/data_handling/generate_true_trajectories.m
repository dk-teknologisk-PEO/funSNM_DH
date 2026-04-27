function true_trajectories = generate_true_trajectories(ground_truth, timestamps, drift_config)
%GENERATE_TRUE_TRAJECTORIES Creates time-varying ground truth trajectories.
%
%   Only the offset supports drift, and only for a single specified house.
%   All other houses get constant trajectories.
%
%   Args:
%       ground_truth (table): Ground truth table with columns
%           house_id, true_offset, true_U.
%       timestamps (datetime array): All timesteps.
%       drift_config (struct): Drift configuration:
%           .type (char): 'none', 'linear', or 'step'
%           .house_index (int): Index of the house to drift (1-based).
%           For 'linear':
%               .offset_drift_per_year (scalar): Drift rate [°C/year].
%           For 'step':
%               .step_time (datetime): Time of step change.
%               .offset_step (scalar): Step magnitude [°C].
%
%   Returns:
%       true_trajectories (cell): Cell array of structs, one per house.
%           Each struct has .offset (1xT) and .U (1xT).

    num_houses = height(ground_truth);
    T = length(timestamps);
    true_trajectories = cell(num_houses, 1);

    % Time in years from start
    t0 = timestamps(1);
    t_years = years(timestamps - t0);

    % Determine which house drifts (-1 means none)
    if strcmp(drift_config.type, 'none') || ~isfield(drift_config, 'house_index')
        drift_index = -1;
    else
        drift_index = drift_config.house_index;
    end

    for i = 1:num_houses
        base_offset = ground_truth.true_offset(i);
        base_U = ground_truth.true_U(i);

        % U-value is always constant
        traj_U = repmat(base_U, 1, T);

        % Default: constant offset
        traj_offset = repmat(base_offset, 1, T);

        % Apply drift only to the specified house
        if i == drift_index
            switch drift_config.type
                case 'linear'
                    traj_offset = base_offset + drift_config.offset_drift_per_year * t_years;

                case 'step'
                    after_step = timestamps >= drift_config.step_time;
                    traj_offset(after_step) = base_offset + drift_config.offset_step;
            end
        end

        true_trajectories{i}.offset = traj_offset;
        true_trajectories{i}.U = traj_U;
    end
end