function meter_data_csac = apply_offset_drift_to_data(meter_data_csac, house_ids, drift_config, timestamps)
%APPLY_OFFSET_DRIFT_TO_DATA Applies time-varying offset drift to meter data.
%
%   The base offset has already been applied during initialization.
%   This function applies the ADDITIONAL drift relative to the base offset.
%
%   Args:
%       meter_data_csac (table): Meter data (already has base offset applied).
%       house_ids (1xN int): House IDs.
%       drift_config (struct): Drift configuration.
%       timestamps (datetime array): All unique timestamps (for t0 reference).
%
%   Returns:
%       meter_data_csac (table): Meter data with drift applied.

    if strcmp(drift_config.type, 'none')
        return;
    end

    t0 = timestamps(1);

    switch drift_config.type
        case 'linear'
            t_years = years(meter_data_csac.timestamp - t0);
            additional_offset = drift_config.offset_drift_per_year * t_years;
            meter_data_csac.T_supply_C = meter_data_csac.T_supply_C - additional_offset;
            fprintf('  Applied linear offset drift: %.3f °C/year to %d houses\n', ...
                drift_config.offset_drift_per_year, numel(house_ids));

        case 'step'
            after_step = meter_data_csac.timestamp >= drift_config.step_time;
            meter_data_csac.T_supply_C(after_step) = ...
                meter_data_csac.T_supply_C(after_step) - drift_config.offset_step;
            fprintf('  Applied step offset change: %.3f °C at %s to %d houses\n', ...
                drift_config.offset_step, string(drift_config.step_time), numel(house_ids));

        otherwise
            error('Unknown drift type: %s', drift_config.type);
    end
end