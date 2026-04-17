function [season_state, actions] = manage_heating_season(season_state, current_date, current_timestep, daily_T_air_max_table, config, csac_id)
%MANAGE_HEATING_SEASON Evaluates the heating season gate and returns actions.
%
%   This function encapsulates the full heating season state machine:
%   - Checks daily max air temperature against a lookback window
%   - Detects season starts and ends
%   - Manages cooldown periods after real seasons
%   - Determines when snapshots should be saved or restored
%
%   Args:
%       season_state (struct): Persistent state of the season gate
%           (created by initialize_season_state).
%       current_date (datetime): The current date (start of day).
%       current_timestep (int): The current timestep index.
%       daily_T_air_max_table (table): Precomputed table with columns
%           'date' and 'T_air_max'.
%       config (struct): Project configuration struct.
%       csac_id (int): CSAC identifier for log messages.
%
%   Returns:
%       season_state (struct): Updated season gate state.
%       actions (struct): Actions to be performed by the caller:
%           .season_is_active (logical): Whether processing should proceed.
%           .do_restore_snapshot (logical): Whether to restore states from snapshot.
%           .gap_days (double): Gap duration for P inflation (only valid if restoring).
%           .do_save_snapshot (logical): Whether to save current states as snapshot.
%           .skip_timestep (logical): Whether to skip processing this timestep.

    %% Unpack config
    gate_cfg = config.project.heating_season_gate;
    T_air_max_threshold = gate_cfg.T_air_max_threshold;
    lookback_days = gate_cfg.lookback_days;
    min_cooldown_days = gate_cfg.min_cooldown_days;
    min_season_for_cooldown = gate_cfg.min_season_days_for_cooldown;

    %% Initialize actions (defaults: do nothing special)
    actions.season_is_active = season_state.active;
    actions.do_restore_snapshot = false;
    actions.gap_days = 0;
    actions.do_save_snapshot = false;
    actions.skip_timestep = false;

    %% Only evaluate once per day
    if ~isnat(season_state.last_gate_check_date) && current_date <= season_state.last_gate_check_date
        % Already checked today — return current state
        if ~season_state.active
            actions.skip_timestep = true;
        end
        return;
    end

    % New day — evaluate the gate
    season_state.last_gate_check_date = current_date;

    [gate_is_active, recent_temps, ~] = check_heating_season_gate(...
        current_date, lookback_days, T_air_max_threshold, daily_T_air_max_table);

    was_active = season_state.active;

    %% === TRANSITION: Inactive -> Active (Season Start) ===
    if gate_is_active && ~was_active
        allow_restart = evaluate_cooldown(season_state, current_date, min_cooldown_days, csac_id);

        if allow_restart
            % Clear cooldown if it expired (helper detected but can't modify state)
            if season_state.cooldown_active && ~isnat(season_state.last_real_season_end_date)
                days_since_real_end = days(current_date - season_state.last_real_season_end_date);
                if days_since_real_end >= min_cooldown_days
                    season_state.cooldown_active = false;
                end
            end

            season_state.active = true;
            season_state.count = season_state.count + 1;
            season_state.start_timestep = current_timestep;
            season_state.start_date = current_date;

            fprintf('\n*** CSAC %d: HEATING SEASON #%d START at %s (timestep %d) ***\n', ...
                csac_id, season_state.count, string(current_date), current_timestep);
            fprintf('    Last %d days T_air_max: [%s] (all < %.0f°C)\n', ...
                lookback_days, ...
                format_temp_list(recent_temps), ...
                T_air_max_threshold);

            if season_state.count > 1
                if ~isnat(season_state.last_active_date)
                    actions.gap_days = days(current_date - season_state.last_active_date);
                else
                    actions.gap_days = 0;
                end
                actions.do_restore_snapshot = true;
                fprintf('    Restarting from snapshot (timestep %d), gap = %d days\n', ...
                    season_state.snapshot_timestep, actions.gap_days);
            else
                fprintf('    First season — using initial P (cold start)\n');
            end
        end
    %% === TRANSITION: Active -> Inactive (Season End) ===
    elseif ~gate_is_active && was_active
        season_state.active = false;
        season_state.end_timestep = current_timestep;

        % Calculate season duration
        if ~isnat(season_state.start_date)
            season_state.last_season_duration_days = days(current_date - season_state.start_date);
        else
            season_state.last_season_duration_days = 0;
        end

        fprintf('\n*** CSAC %d: HEATING SEASON #%d END at %s (timestep %d, ran %d days) ***\n', ...
            csac_id, season_state.count, string(current_date), current_timestep, ...
            season_state.last_season_duration_days);
        fprintf('    T_air_max exceeded %.0f°C in lookback window: [%s]\n', ...
            T_air_max_threshold, ...
            format_temp_list(recent_temps));

        % Only activate cooldown and flag snapshot save for real seasons
        if season_state.last_season_duration_days >= min_season_for_cooldown
            actions.do_save_snapshot = true;
            season_state.last_real_season_end_date = current_date;
            season_state.cooldown_active = true;
            season_state.snapshot_timestep = current_timestep;
            fprintf('    Snapshot flagged for save (real season, cooldown activated)\n');
        else
            fprintf('    Snapshot NOT saved (season only ran %d days < %d minimum)\n', ...
                season_state.last_season_duration_days, min_season_for_cooldown);
            if season_state.cooldown_active
                fprintf('    Cooldown still active from real season ending %s\n', ...
                    string(season_state.last_real_season_end_date));
            end
        end
    end

    %% Update gate statistics
    if season_state.active
        season_state.active_days = season_state.active_days + 1;
    else
        season_state.inactive_days = season_state.inactive_days + 1;
    end

    %% Set output actions
    actions.season_is_active = season_state.active;
    if ~season_state.active
        actions.skip_timestep = true;
    end
end


%% ===== LOCAL HELPER FUNCTIONS =====

function allow = evaluate_cooldown(season_state, current_date, min_cooldown_days, csac_id)
%EVALUATE_COOLDOWN Checks whether cooldown permits a season restart.
    allow = true;

    if ~season_state.cooldown_active || isnat(season_state.last_real_season_end_date)
        return;
    end

    days_since_real_end = days(current_date - season_state.last_real_season_end_date);
    current_month_num = month(current_date);
    is_winter_month = (current_month_num >= 10) || (current_month_num <= 3);

    if days_since_real_end >= min_cooldown_days
        % Full cooldown elapsed — will be cleared by caller via season_state update
        % (We don't modify season_state here since it's not passed by reference;
        %  the caller clears cooldown_active when allow==true and cooldown was active.)
        fprintf('    CSAC %d: Cooldown expired (%d days since %s)\n', ...
            csac_id, round(days_since_real_end), string(season_state.last_real_season_end_date));
    elseif is_winter_month
        % Winter month — allow restart but keep cooldown active
        fprintf('    CSAC %d: Cooldown bypassed (winter month %d), cooldown remains active\n', ...
            csac_id, current_month_num);
    else
        % Spring/summer and cooldown hasn't elapsed — block
        allow = false;
        if mod(season_state.inactive_days, 30) == 1
            fprintf('    CSAC %d: Cooldown active — %d/%d days since real season end at %s\n', ...
                csac_id, round(days_since_real_end), min_cooldown_days, ...
                string(season_state.last_real_season_end_date));
        end
    end
end


function s = format_temp_list(temps)
%FORMAT_TEMP_LIST Formats a vector of temperatures as a comma-separated string.
    s = strjoin(arrayfun(@(x) sprintf('%.1f', x), temps(isfinite(temps)), 'UniformOutput', false), ', ');
end