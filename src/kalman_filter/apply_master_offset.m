function ukf_states = apply_master_offset(ukf_states, master_offset, num_active_houses, config)
%APPLY_MASTER_OFFSET Applies a damped, clamped master offset to all house offsets.
%
%   The master offset is a shared bias detected during T_main fitting.
%   It is applied to all houses to correct systematic errors.
%
%   Args:
%       ukf_states (cell): UKF state structs for all houses.
%       master_offset (scalar): Raw master offset from T_main fitting.
%       num_active_houses (int): Number of houses with sufficient flow.
%       config (struct): Project configuration.
%
%   Returns:
%       ukf_states (cell): UKF states with adjusted offsets.

    gamma = config.project.master_offset.gamma;
    deadzone = config.project.master_offset.deadzone;
    min_active = config.project.initialization.min_active_houses;
    max_master_offset = 0.5;

    if isnan(master_offset) || abs(master_offset) <= deadzone
        return;
    end

    if num_active_houses < min_active
        return;
    end

    clamped = max(-max_master_offset, min(max_master_offset, master_offset));
    damped = gamma * clamped;

    for i = 1:numel(ukf_states)
        ukf_states{i}.x(1) = ukf_states{i}.x(1) - damped;
    end
end