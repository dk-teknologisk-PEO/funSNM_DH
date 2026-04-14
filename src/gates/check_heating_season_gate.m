function [is_active, recent_max_temps, n_valid] = check_heating_season_gate(...
        current_date, lookback_days, threshold, daily_T_air_max_tbl)
%CHECK_HEATING_SEASON_GATE Check if heating season is active based on temperature.
%   Returns true if ALL of the last N days had T_air_max < threshold.
    
    recent_max_temps = nan(lookback_days, 1);
    for d = 0:(lookback_days - 1)
        check_date = current_date - days(d);
        row_idx = find(daily_T_air_max_tbl.date == check_date, 1);
        if ~isempty(row_idx)
            recent_max_temps(d + 1) = daily_T_air_max_tbl.T_air_max(row_idx);
        end
    end
    
    n_valid = sum(isfinite(recent_max_temps));
    
    if n_valid < lookback_days
        is_active = false;
        return;
    end
    
    % ALL days in the lookback window must be below threshold
    is_active = all(recent_max_temps(isfinite(recent_max_temps)) < threshold);
end