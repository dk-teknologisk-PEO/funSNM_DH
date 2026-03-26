% save_logger_to_csv.m

function save_logger_to_csv(logger, output_folder, file_prefix)
%SAVE_LOGGER_TO_CSV Converts the logger structure to a tidy CSV file.
%
%   Args:
%       logger (struct): The logger structure from the simulation.
%       output_folder (string): The path to the folder to save the CSV.
%       file_prefix (string): A prefix for the filename (e.g., 'ukf' or 'pf').

    num_houses = length(logger.house_ids);
    num_timestamps = length(logger.timestamps);
    
    % Pre-allocate a cell array to hold all data before converting to table
    % This is much faster than building the table row-by-row
    all_data = cell(num_timestamps * num_houses, 10);
    
    rowIndex = 1;
    for t = 1:num_timestamps
        timestamp_val = logger.timestamps(t);
        % Skip if timestamp is not valid (NaT)
        if isnat(timestamp_val)
            continue;
        end
        
        for h = 1:num_houses
            house_id_val = logger.house_ids(h);
            
            all_data{rowIndex, 1} = timestamp_val;
            all_data{rowIndex, 2} = house_id_val;
            
            % State Estimates
            all_data{rowIndex, 3} = logger.state_estimates(1, h, t); % Offset
            all_data{rowIndex, 4} = logger.state_estimates(2, h, t); % U-Value
            
            % Uncertainty (Standard Deviation)
            all_data{rowIndex, 5} = sqrt(logger.covariance_posterior(1, h, t)); % StdDev Offset
            all_data{rowIndex, 6} = sqrt(logger.covariance_posterior(2, h, t)); % StdDev U-Value
            
            % Diagnostics
            all_data{rowIndex, 7} = logger.innovations(h, t);
            all_data{rowIndex, 8} = logger.innovation_variance(h, t);
            all_data{rowIndex, 9} = logger.nis(h, t);
            
            % Kalman Gain (just the first component for simplicity)
            all_data{rowIndex, 10} = logger.kalman_gain(1, h, t);
            
            rowIndex = rowIndex + 1;
        end
    end
    
    % Trim any unused rows if there were NaTs
    all_data = all_data(1:rowIndex-1, :);
    
    % Create the table with descriptive variable names
    col_names = {'Timestamp', 'HouseID', 'Offset_Estimate_C', 'UValue_Estimate_W_m_K', ...
                 'Offset_StdDev_C', 'UValue_StdDev_W_m_K', 'Innovation_C', ...
                 'Innovation_Variance', 'NIS', 'KalmanGain_Offset'};
    
    output_table = cell2table(all_data, 'VariableNames', col_names);
    
    % Sort by time and house for consistency
    output_table = sortrows(output_table, {'Timestamp', 'HouseID'});
    
    % Write to CSV
    filename = fullfile(output_folder, strcat(file_prefix, '_filter_outputs.csv'));
    writetable(output_table, filename);
    
    fprintf('Saved filter log to %s\n', filename);

end