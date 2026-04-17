function print_csac_summary(csac_id, gate_accept_count, gate_reject_count, season_state)
%PRINT_CSAC_SUMMARY Prints final statistics for a CSAC after processing.
%
%   Args:
%       csac_id (int): CSAC identifier.
%       gate_accept_count (int): Total accepted UKF updates.
%       gate_reject_count (int): Total rejected UKF updates.
%       season_state (struct): Final heating season state.

    fprintf('\n=== CSAC %d Final Statistics ===\n', csac_id);
    fprintf('UKF: %d accepted, %d rejected (%.1f%% rejection rate)\n', ...
        gate_accept_count, gate_reject_count, ...
        100 * gate_reject_count / max(1, gate_accept_count + gate_reject_count));
    fprintf('Season gate: %d active days, %d inactive days (%.1f%% active)\n', ...
        season_state.active_days, season_state.inactive_days, ...
        100 * season_state.active_days / max(1, season_state.active_days + season_state.inactive_days));
    fprintf('Heating seasons started: %d\n', season_state.count);
    if ~isnan(season_state.end_timestep)
        fprintf('Last season ended at timestep %d\n', season_state.end_timestep);
    end
    if season_state.cooldown_active
        fprintf('Cooldown active from real season ending %s\n', ...
            string(season_state.last_real_season_end_date));
    end
    fprintf('=====================================\n\n');
end