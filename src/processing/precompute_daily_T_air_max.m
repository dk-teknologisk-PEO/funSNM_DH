function daily_T_air_max_table = precompute_daily_T_air_max(T_air_C)
%PRECOMPUTE_DAILY_T_AIR_MAX Build daily max air temperature lookup table.

    fprintf('Pre-computing daily T_air_max lookup table...\n');
    
    T_air_C.date = dateshift(T_air_C.time, 'start', 'day');
    [air_date_groups, air_unique_dates] = findgroups(T_air_C.date);
    daily_T_air_max_values = splitapply(@max, T_air_C.values, air_date_groups);
    
    daily_T_air_max_table = table(air_unique_dates, daily_T_air_max_values, ...
        'VariableNames', {'date', 'T_air_max'});
    daily_T_air_max_table = sortrows(daily_T_air_max_table, 'date');
    
    fprintf('  Built daily T_air_max table: %d days (%.1f to %.1f °C)\n', ...
        height(daily_T_air_max_table), ...
        min(daily_T_air_max_table.T_air_max), ...
        max(daily_T_air_max_table.T_air_max));
end