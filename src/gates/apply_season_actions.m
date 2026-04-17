function [ukf_states, ukf_states_snapshot] = apply_season_actions(actions, ukf_states, ukf_states_snapshot, P_base, config)
%APPLY_SEASON_ACTIONS Executes snapshot restore/save actions from the season gate.
%
%   This function modifies UKF states based on the actions struct returned
%   by manage_heating_season. It handles:
%   - Restoring states from snapshot with gap-proportional P inflation
%   - Saving current states as a new snapshot
%
%   Args:
%       actions (struct): Actions struct from manage_heating_season.
%       ukf_states (cell array): Current UKF state structs for each house.
%       ukf_states_snapshot (cell array): Snapshot UKF state structs.
%       P_base (2x2 matrix): Initial covariance (cold start upper bound).
%       config (struct): Project configuration.
%
%   Returns:
%       ukf_states (cell array): Possibly modified UKF states.
%       ukf_states_snapshot (cell array): Possibly updated snapshot.

    gate_cfg = config.project.heating_season_gate;
    gap_P_growth_offset = gate_cfg.gap_P_growth_per_day_offset;
    gap_P_growth_U = gate_cfg.gap_P_growth_per_day_U;
    num_houses = numel(ukf_states);

    if actions.do_restore_snapshot
        gap_days = actions.gap_days;
        for i = 1:num_houses
            % Restore point estimates from snapshot
            ukf_states{i}.x = ukf_states_snapshot{i}.x;

            % Grow P proportionally to gap duration
            P_snap = ukf_states_snapshot{i}.P;
            P_grown = P_snap;
            P_grown(1,1) = P_snap(1,1) + (gap_P_growth_offset^2) * gap_days;
            P_grown(2,2) = P_snap(2,2) + (gap_P_growth_U^2) * gap_days;

            % Cap at P_base (never more uncertain than cold start)
            P_grown(1,1) = min(P_grown(1,1), P_base(1,1));
            P_grown(2,2) = min(P_grown(2,2), P_base(2,2));

            ukf_states{i}.P = P_grown;

            if i == 1
                fprintf('    House 1: P_snap=[%.4f, %.6f], P_grown=[%.4f, %.6f], P_base=[%.4f, %.6f]\n', ...
                    sqrt(P_snap(1,1)), sqrt(P_snap(2,2)), ...
                    sqrt(P_grown(1,1)), sqrt(P_grown(2,2)), ...
                    sqrt(P_base(1,1)), sqrt(P_base(2,2)));
            end
        end
    end

    if actions.do_save_snapshot
        for i = 1:num_houses
            ukf_states_snapshot{i}.x = ukf_states{i}.x;
            ukf_states_snapshot{i}.P = ukf_states{i}.P;
        end
    end
end