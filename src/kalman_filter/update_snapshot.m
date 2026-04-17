function [ukf_states_snapshot, snapshot_timestep] = update_snapshot(season_state, current_date, t, ukf_states, ukf_states_snapshot, config)
%UPDATE_SNAPSHOT Saves a snapshot of UKF states if the season is established.
%
%   Only saves when the current season has run longer than the minimum
%   season duration required for cooldown activation.
%
%   Args:
%       season_state (struct): Current heating season state.
%       current_date (datetime): Current date.
%       t (int): Current timestep index.
%       ukf_states (cell): Current UKF states.
%       ukf_states_snapshot (cell): Current snapshot states.
%       config (struct): Project configuration.
%
%   Returns:
%       ukf_states_snapshot (cell): Possibly updated snapshot.
%       snapshot_timestep (double): Timestep of snapshot (NaN or updated t).

    snapshot_timestep = season_state.snapshot_timestep;
    min_season_for_cooldown = config.project.heating_season_gate.min_season_days_for_cooldown;

    if ~isnat(season_state.start_date)
        current_season_days = days(current_date - season_state.start_date);
    else
        current_season_days = 0;
    end

    if current_season_days >= min_season_for_cooldown
        snapshot_timestep = t;
        for i = 1:numel(ukf_states)
            ukf_states_snapshot{i}.x = ukf_states{i}.x;
            ukf_states_snapshot{i}.P = ukf_states{i}.P;
        end
    end
end