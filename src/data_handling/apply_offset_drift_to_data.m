function meter_data_csac = apply_offset_drift_to_data(meter_data_csac, house_ids, drift_config, timestamps)
%APPLY_OFFSET_DRIFT_TO_DATA Applies time-varying offset drift to meter data.
%
%   The base offset has already been applied during initialization.
%   This function applies ADDITIONAL drift to a single specified house.
%
%   Args:
%       meter_data_csac (table): Meter data (already has base offset applied).
%       house_ids (1xN int): House IDs in this CSAC.
%       drift_config (struct): Drift configuration:
%           .type (char): 'none', 'linear', or 'step'
%           .house_index (int): Index into house_ids of the house to drift.
%           For 'linear':
%               .offset_drift_per_year (scalar): Drift rate [°C/year].
%           For 'step':
%               .step_time (datetime): Time of step change.
%               .offset_step (scalar): Step magnitude [°C].
%
%   Returns:
%       meter_data_csac (table): Meter data with drift applied to one house.

    if strcmp(drift_config.type, 'none')
        return;
    end

    t0 = timestamps(1);
    target_house_id = house_ids(drift_config.house_index);
    mask = meter_data_csac.house_id == target_house_id;

    switch drift_config.type
        case 'linear'
            t_years = years(meter_data_csac.timestamp(mask) - t0);
            additional_offset = drift_config.offset_drift_per_year * t_years;
            meter_data_csac.T_supply_C(mask) = meter_data_csac.T_supply_C(mask) - additional_offset;
            fprintf('  Applied linear offset drift: %.3f °C/year to house %d (index %d)\n', ...
                drift_config.offset_drift_per_year, target_house_id, drift_config.house_index);

        case 'step'
            step_mask = mask & (meter_data_csac.timestamp >= drift_config.step_time);
            meter_data_csac.T_supply_C(step_mask) = ...
                meter_data_csac.T_supply_C(step_mask) - drift_config.offset_step;
            fprintf('  Applied step offset change: %.3f °C at %s to house %d (index %d)\n', ...
                drift_config.offset_step, string(drift_config.step_time), ...
                target_house_id, drift_config.house_index);

        otherwise
            error('Unknown drift type: %s', drift_config.type);
    end
end