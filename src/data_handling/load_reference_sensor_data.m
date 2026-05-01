function ref_data = load_reference_sensor_data(network_data, sensor_config, network_id)
%LOAD_REFERENCE_SENSOR_DATA Extracts reference sensor temperatures and adds noise.
%
%   Args:
%       network_data (table): Network temperature data with columns like
%           Plant_s, J_Main_0_s, J_Main_1_s, etc.
%       sensor_config (struct): Sensor configuration:
%           .type (char): 'none', 'plant', or 'csac_junction'
%           .csac_id (int): CSAC ID for junction sensor (ignored for 'plant')
%           .noise_std (scalar): Measurement noise std [°C]
%           .seed (int): Random seed for reproducible noise
%       network_id (int): Network identifier for seeding.
%
%   Returns:
%       ref_data (struct):
%           .timestamps (datetime array): Timestamps.
%           .temperatures (double array): Noisy reference temperatures.
%           .location (char): Description of sensor location.
%           .column_name (char): Original column name used.
%           .valid (logical): Whether sensor data was found.

    ref_data.valid = false;
    ref_data.timestamps = network_data.timestamp;
    ref_data.temperatures = nan(height(network_data), 1);
    ref_data.location = 'none';
    ref_data.column_name = '';

    if strcmp(sensor_config.type, 'none')
        return;
    end

    %% Determine which column to use
    switch sensor_config.type
        case 'plant'
            col_name = 'Plant_s';
            location = 'Network inlet';

        case 'csac_junction'
            col_name = sprintf('J_Main_%d_s', sensor_config.csac_id);
            location = sprintf('CSAC %d junction on main pipe', sensor_config.csac_id);

        otherwise
            warning('Unknown sensor type: %s', sensor_config.type);
            return;
    end

    %% Check column exists
    if ~ismember(col_name, network_data.Properties.VariableNames)
        warning('Column %s not found in network data', col_name);
        return;
    end

    %% Extract temperatures and add noise
    true_temps = network_data.(col_name);

    % Reproducible noise
    rng(network_id * 10000 + 5e6 + sensor_config.csac_id, 'twister');
    noise = sensor_config.noise_std * randn(size(true_temps));
    rng('shuffle');

    ref_data.temperatures = true_temps + noise;
    ref_data.location = location;
    ref_data.column_name = col_name;
    ref_data.valid = true;

    fprintf('  Reference sensor loaded: %s (%s), noise=%.2f°C\n', ...
        location, col_name, sensor_config.noise_std);
end