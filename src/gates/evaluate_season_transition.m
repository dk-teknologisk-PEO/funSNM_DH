function csac_state = evaluate_season_transition(csac_state, t, current_date, ...
        daily_T_air_max_table, params)
%EVALUATE_SEASON_TRANSITION Check for heating season start/end transitions.
%   Also handles gradual P growth during dormant periods and cumulative
%   winter duration tracking for cooldown.

    %% Gradual uncertainty growth when dormant
    % Add process noise daily regardless of season state, so that
    % uncertainty grows naturally during gaps rather than in a step function
    if ~csac_state.season_active
        if isnat(csac_state.last_dormant_P_growth_date) || current_date > csac_state.last_dormant_P_growth_date
            csac_state.last_dormant_P_growth_date = current_date;
            
            % Add one day's worth of process noise (scaled up since updates
            % normally happen ~4x/day, so one dormant day ≈ 4 missed steps)
            steps_per_day = 4;
            % Determine P ceiling based on whether we've had a successful season
            if csac_state.cumulative_winter_days >= params.season_gate_min_season_for_cooldown
                P_ceil = params.P_max_restart;
            else
                P_ceil = params.P_base;  % First winter — allow full cold-start range
            end
            
            for i = 1:csac_state.num_houses
                P_grown = csac_state.ukf_states{i}.P + steps_per_day * params.Q_base;
                P_grown(1,1) = min(P_grown(1,1), P_ceil(1,1));
                P_grown(2,2) = min(P_grown(2,2), P_ceil(2,2));
                csac_state.ukf_states{i}.P = P_grown;
                
                P_grown_pf = csac_state.pf_states{i}.P + steps_per_day * params.Q_base;
                P_grown_pf(1,1) = min(P_grown_pf(1,1), P_ceil(1,1));
                P_grown_pf(2,2) = min(P_grown_pf(2,2), P_ceil(2,2));
                csac_state.pf_states{i}.P = P_grown_pf;
            end
        end
    end

    %% Daily gate check
    if ~(isnat(csac_state.last_gate_check_date) || current_date > csac_state.last_gate_check_date)
        return;
    end
    
    csac_state.last_gate_check_date = current_date;
    
    [gate_is_active, recent_temps, ~] = check_heating_season_gate(...
        current_date, params.season_gate_lookback_days, ...
        params.season_gate_T_air_max_threshold, daily_T_air_max_table);
    
    csac_state.was_season_active = csac_state.season_active;
    
    if gate_is_active && ~csac_state.was_season_active
        %% === Potential season start ===
        allow_restart = true;
        
        % Check cooldown using CUMULATIVE winter duration
        if ~isnat(csac_state.last_season_end_date) && ...
                csac_state.cumulative_winter_days >= params.season_gate_min_season_for_cooldown
            days_since_end = days(current_date - csac_state.last_season_end_date);
            current_month = month(current_date);
            is_winter_month = (current_month >= 10) || (current_month <= 3);
            
            if ~is_winter_month && days_since_end < params.season_gate_min_cooldown_days
                allow_restart = false;
                if mod(csac_state.season_gate_inactive_days, 30) == 1
                    fprintf('    Cooldown active — %d/%d days since season end (cumulative winter: %d days)\n', ...
                        round(days_since_end), params.season_gate_min_cooldown_days, ...
                        csac_state.cumulative_winter_days);
                end
            end
        end
        
        if allow_restart
            csac_state.season_active = true;
            csac_state.season_count = csac_state.season_count + 1;
            csac_state.season_start_timestep = t;
            csac_state.season_start_date = current_date;
            
            % Reset cumulative counter at start of a new winter (Oct-Nov start)
            current_month = month(current_date);
            if current_month >= 10 && current_month <= 11
                csac_state.cumulative_winter_days = 0;
                fprintf('\n*** HEATING SEASON #%d START at %s (timestep %d) — new winter ***\n', ...
                    csac_state.season_count, string(current_date), t);
            else
                fprintf('\n*** HEATING SEASON #%d START at %s (timestep %d) — continuing winter ***\n', ...
                    csac_state.season_count, string(current_date), t);
            end
            
            fprintf('    Last %d days T_air_max: [%s] (all < %.0f°C)\n', ...
                params.season_gate_lookback_days, ...
                strjoin(arrayfun(@(x) sprintf('%.1f', x), ...
                recent_temps(isfinite(recent_temps)), 'UniformOutput', false), ', '), ...
                params.season_gate_T_air_max_threshold);
            fprintf('    Cumulative winter days so far: %d\n', csac_state.cumulative_winter_days);
            
            if csac_state.season_count == 1
                fprintf('    First season — using initial P (cold start)\n');
            else
                % No big P inflation — uncertainty has been growing gradually
                % Just restore point estimates from snapshot
                fprintf('    Restarting from snapshot (timestep %d)\n', csac_state.snapshot_timestep);
                fprintf('    P has grown gradually during dormancy — no step inflation\n');
                for i = 1:csac_state.num_houses
                    csac_state.ukf_states{i}.x = csac_state.ukf_states_snapshot{i}.x;
                    csac_state.pf_states{i}.x = csac_state.pf_states_snapshot{i}.x;
                    % P is already at the right level from gradual growth
                    % Don't overwrite it — just keep current P
                    
                    if i == 1
                        fprintf('    House %d: current P=[%.4f, %.6f]\n', ...
                            csac_state.house_ids(i), ...
                            sqrt(csac_state.ukf_states{i}.P(1,1)), ...
                            sqrt(csac_state.ukf_states{i}.P(2,2)));
                    end
                end
            end
        end
        
    elseif ~gate_is_active && csac_state.was_season_active
        %% === Season end ===
        csac_state.season_active = false;
        csac_state.season_end_timestep = t;
        csac_state.last_season_end_date = current_date;
        
        % Calculate this segment's duration and add to cumulative
        if ~isnat(csac_state.season_start_date)
            segment_days = days(current_date - csac_state.season_start_date);
        else
            segment_days = 0;
        end
        csac_state.last_season_duration_days = segment_days;
        csac_state.cumulative_winter_days = csac_state.cumulative_winter_days + segment_days;
        
        fprintf('\n*** HEATING SEASON #%d END at %s (timestep %d, ran %d days, cumulative: %d days) ***\n', ...
            csac_state.season_count, string(current_date), t, segment_days, ...
            csac_state.cumulative_winter_days);
        
        % Save snapshot if cumulative winter is long enough
        if csac_state.cumulative_winter_days >= params.season_gate_min_season_for_cooldown
            for i = 1:csac_state.num_houses
                csac_state.ukf_states_snapshot{i}.x = csac_state.ukf_states{i}.x;
                csac_state.ukf_states_snapshot{i}.P = csac_state.ukf_states{i}.P;
                csac_state.pf_states_snapshot{i}.x = csac_state.pf_states{i}.x;
                csac_state.pf_states_snapshot{i}.P = csac_state.pf_states{i}.P;
            end
            csac_state.snapshot_timestep = t;
            fprintf('    Snapshot saved (cumulative winter >= %d days)\n', ...
                params.season_gate_min_season_for_cooldown);
        else
            fprintf('    Snapshot NOT saved (cumulative winter only %d days < %d minimum)\n', ...
                csac_state.cumulative_winter_days, params.season_gate_min_season_for_cooldown);
        end
    end
    
    %% Track gate statistics
    if csac_state.season_active
        csac_state.season_gate_active_days = csac_state.season_gate_active_days + 1;
    else
        csac_state.season_gate_inactive_days = csac_state.season_gate_inactive_days + 1;
    end
end