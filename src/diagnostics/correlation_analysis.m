% correlation_analysis.m
% Analyzes correlation between heating demand stability and weather variables
% to determine optimal heating season gating criteria.

clear all
close all
clc

%% ============================================================
%% 1. CONFIGURATION
%% ============================================================

% Paths — adjust these to your setup
utility_data_path = 'C:\Users\peo\Documents\GitHub\DH simulation\data\processed_data_folder';      % folder containing utility CSV files
weather_data_path = 'C:\Users\peo\Documents\GitHub\DH simulation\data\weather';      % folder containing weather CSV files

% Weather file names (adjust to match your files)
air_temp_file    = fullfile(weather_data_path, 'mean_temp.csv');
soil_temp_file   = fullfile(weather_data_path, 'temp_soil_30.csv');
radiation_file   = fullfile(weather_data_path, 'mean_radiation.csv');
cloud_cover_file = fullfile(weather_data_path, 'mean_cloud_cover.csv');
wind_speed_file  = fullfile(weather_data_path, 'mean_wind_speed.csv');

% Time range
t_start = datetime(2018, 1, 1);
t_end   = datetime(2020, 12, 31, 23, 0, 0);

% Minimum flow to consider a house "active" in an hour
min_hourly_flow_m3 = 0.005;  % 5 liters/hour

fprintf('=== District Heating Correlation Analysis ===\n');
fprintf('Period: %s to %s\n\n', string(t_start), string(t_end));


%% ============================================================
%% 3. LOAD WEATHER DATA
%% ============================================================

fprintf('Loading weather data...\n');

% --- Air temperature ---
if isfile(air_temp_file)
    air_hourly = load_weather_file(air_temp_file, 'T_air_C');
else
    error('Air temperature file not found: %s', air_temp_file);
end

% --- Soil temperature ---
if isfile(soil_temp_file)
    soil_hourly = load_weather_file(soil_temp_file, 'T_soil_C');
else
    warning('Soil temperature file not found, skipping.');
    soil_hourly = table();
end

% --- Solar radiation ---
if isfile(radiation_file)
    rad_hourly = load_weather_file(radiation_file, 'radiation_W_m2');
else
    warning('Radiation file not found, skipping.');
    rad_hourly = table();
end

% --- Cloud cover ---
if isfile(cloud_cover_file)
    cloud_hourly = load_weather_file(cloud_cover_file, 'cloud_cover');
else
    warning('Cloud cover file not found, skipping.');
    cloud_hourly = table();
end

% --- Wind speed ---
if isfile(wind_speed_file)
    wind_hourly = load_weather_file(wind_speed_file, 'wind_speed_m_s');
else
    warning('Wind speed file not found, skipping.');
    wind_hourly = table();
end

%% ============================================================
%% 4. LOAD AND PROCESS UTILITY DATA
%% ============================================================

fprintf('\nLoading utility data...\n');

% Find all CSV files in the utility data folder
utility_files = dir(fullfile(utility_data_path, '*.csv'));
fprintf('  Found %d utility data files\n', length(utility_files));

% We'll accumulate hourly power per house, then aggregate daily
all_hourly_power = [];

num_files_processed = 0;
num_houses_total = 0;
num_houses_to_analyze = 500;

num_houses_to_analyze = min(num_houses_to_analyze, length(utility_files));

utility_files_idx = randperm(length(utility_files));

for f = 1:num_houses_to_analyze
    idx = utility_files_idx(f);
    filepath = fullfile(utility_files(idx).folder, utility_files(idx).name);
    
    try
        raw = readtable(filepath, 'TextType', 'string');
    catch
        warning('Could not read file: %s', utility_files(idx).name);
        continue;
    end
    
    % Check required columns exist
    required_cols = {'customer_id', 'time_rounded', 'energy_heat_kwh', 'volume_flow_m3'};
    if ~all(ismember(required_cols, raw.Properties.VariableNames))
        warning('Missing columns in %s, skipping.', utility_files(idx).name);
        continue;
    end
    
    % Ensure numeric types
    if ~isnumeric(raw.energy_heat_kwh)
        raw.energy_heat_kwh = str2double(raw.energy_heat_kwh);
    end
    if ~isnumeric(raw.volume_flow_m3)
        raw.volume_flow_m3 = str2double(raw.volume_flow_m3);
    end
    if ~isnumeric(raw.customer_id)
        raw.customer_id = str2double(raw.customer_id);
    end
    
    % Parse timestamps
    try
        raw.time = datetime(raw.time_rounded, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    catch
        try
            raw.time = datetime(raw.time_rounded, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ssXXX');
        catch
            raw.time = datetime(raw.time_rounded);
        end
    end
    if ~isempty(raw.time.TimeZone)
        raw.time.TimeZone = '';
    end
    
    % Filter to time range
    raw = raw(raw.time >= t_start & raw.time <= t_end, :);
    if height(raw) < 48
        continue;
    end
    
    % Get unique houses in this file
    house_ids = unique(raw.customer_id);
    
    for h = 1:length(house_ids)
        hid = house_ids(h);
        hdata = raw(raw.customer_id == hid, :);
        hdata = sortrows(hdata, 'time');
        
        if height(hdata) < 48
            continue;
        end
        
        % Compute hourly differences (cumulative meter readings)
        delta_energy = diff(hdata.energy_heat_kwh);
        delta_volume = diff(hdata.volume_flow_m3);
        delta_time_h = hours(diff(hdata.time));
        
        % Filter out bad intervals
        valid = (delta_time_h > 0.5) & (delta_time_h < 1.5) & ...
                (delta_energy >= 0) & (delta_volume >= 0) & ...
                (delta_energy < 50) & (delta_volume < 5);
        
        n_valid = sum(valid);
        if n_valid < 24
            continue;
        end
        
        % Compute hourly power (kW) and flow rate
        hourly_time = hdata.time(2:end);
        hourly_power_kW = delta_energy ./ delta_time_h;
        hourly_flow_m3_h = delta_volume ./ delta_time_h;
        
        % Apply validity mask
        hourly_time = hourly_time(valid);
        hourly_power_kW = hourly_power_kW(valid);
        hourly_flow_m3_h = hourly_flow_m3_h(valid);
        
        % Store
        n = length(hourly_time);
        house_block = table(hourly_time, ...
            repmat(hid, n, 1), ...
            hourly_power_kW, ...
            hourly_flow_m3_h, ...
            'VariableNames', {'time', 'house_id', 'power_kW', 'flow_m3_h'});
        
        if isempty(all_hourly_power)
            all_hourly_power = house_block;
        else
            all_hourly_power = [all_hourly_power; house_block];
        end
        
        num_houses_total = num_houses_total + 1;
    end
    
    num_files_processed = num_files_processed + 1;
    if mod(num_files_processed, 50) == 0
        fprintf('  Processed %d/%d files, %d houses so far...\n', ...
            num_files_processed, length(utility_files), num_houses_to_analyze);
    end
end

fprintf('  Total: %d houses from %d files, %d hourly records\n', ...
    num_houses_total, num_files_processed, height(all_hourly_power));

if isempty(all_hourly_power)
    error('No valid utility data loaded. Check file paths and column names.');
end

%% ============================================================
%% 5. AGGREGATE TO DAILY METRICS
%% ============================================================

fprintf('\nAggregating to daily metrics...\n');

% Add date column
all_hourly_power.date = dateshift(all_hourly_power.time, 'start', 'day');

% --- Per-house daily aggregation ---
[groups, house_ids_g, dates_g] = findgroups(all_hourly_power.house_id, all_hourly_power.date);

daily_energy   = splitapply(@sum,   all_hourly_power.power_kW,   groups);
daily_mean_pwr = splitapply(@mean,  all_hourly_power.power_kW,   groups);
daily_std_pwr  = splitapply(@std,   all_hourly_power.power_kW,   groups);
daily_n_hours  = splitapply(@numel, all_hourly_power.power_kW,   groups);
daily_active_h = splitapply(@(x) sum(x > 0.1), all_hourly_power.power_kW, groups);
daily_mean_flw = splitapply(@mean,  all_hourly_power.flow_m3_h,  groups);

daily_house = table(house_ids_g, dates_g, daily_energy, daily_mean_pwr, ...
    daily_std_pwr, daily_n_hours, daily_active_h, daily_mean_flw, ...
    'VariableNames', {'house_id', 'date', 'daily_energy_kWh', 'mean_power_kW', ...
                      'std_power_kW', 'n_hours', 'active_hours', 'mean_flow_m3_h'});

% --- Network-wide daily aggregation ---
[date_groups, unique_dates] = findgroups(daily_house.date);

daily_network = table(unique_dates, ...
    splitapply(@sum,   daily_house.daily_energy_kWh, date_groups), ...
    splitapply(@mean,  daily_house.mean_power_kW,    date_groups), ...
    splitapply(@mean,  daily_house.std_power_kW,     date_groups), ...
    splitapply(@numel, daily_house.house_id,         date_groups), ...
    splitapply(@mean,  daily_house.active_hours,     date_groups), ...
    splitapply(@mean,  daily_house.mean_flow_m3_h,   date_groups), ...
    'VariableNames', {'date', 'total_energy_kWh', 'mean_power_kW', 'mean_std_power_kW', ...
                      'n_houses', 'mean_active_hours', 'mean_flow_m3_h'});

% Coefficient of variation of power across houses (daily)
daily_network.cv_power = daily_network.mean_std_power_kW ./ max(daily_network.mean_power_kW, 0.01);

fprintf('  %d unique dates, %d daily-house records\n', height(daily_network), height(daily_house));

%% ============================================================
%% 6. AGGREGATE WEATHER TO DAILY
%% ============================================================

fprintf('\nAggregating weather to daily...\n');

% --- Air temperature: daily mean, max, min ---
air_hourly.date = dateshift(air_hourly.time, 'start', 'day');
[ag, ud] = findgroups(air_hourly.date);

daily_weather = table(ud, ...
    splitapply(@mean,                    air_hourly.T_air_C, ag), ...
    splitapply(@max,                     air_hourly.T_air_C, ag), ...
    splitapply(@min,                     air_hourly.T_air_C, ag), ...
    splitapply(@(x) max(x) - min(x),    air_hourly.T_air_C, ag), ...
    'VariableNames', {'date', 'T_air_mean', 'T_air_max', 'T_air_min', 'T_air_range'});

% --- Soil temperature ---
if ~isempty(soil_hourly)
    soil_hourly.date = dateshift(soil_hourly.time, 'start', 'day');
    [ag, ud] = findgroups(soil_hourly.date);
    soil_daily = table(ud, ...
        splitapply(@mean, soil_hourly.T_soil_C, ag), ...
        'VariableNames', {'date', 'T_soil_mean'});
    daily_weather = outerjoin(daily_weather, soil_daily, 'Keys', 'date', 'MergeKeys', true);
end

% --- Radiation ---
if ~isempty(rad_hourly)
    rad_hourly.date = dateshift(rad_hourly.time, 'start', 'day');
    [ag, ud] = findgroups(rad_hourly.date);
    rad_daily = table(ud, ...
        splitapply(@mean, rad_hourly.radiation_W_m2, ag), ...
        splitapply(@max,  rad_hourly.radiation_W_m2, ag), ...
        splitapply(@sum,  rad_hourly.radiation_W_m2, ag), ...
        'VariableNames', {'date', 'radiation_mean', 'radiation_max', 'radiation_total'});
    daily_weather = outerjoin(daily_weather, rad_daily, 'Keys', 'date', 'MergeKeys', true);
end

% --- Cloud cover ---
if ~isempty(cloud_hourly)
    cloud_hourly.date = dateshift(cloud_hourly.time, 'start', 'day');
    [ag, ud] = findgroups(cloud_hourly.date);
    cloud_daily = table(ud, ...
        splitapply(@mean, cloud_hourly.cloud_cover, ag), ...
        'VariableNames', {'date', 'cloud_cover_mean'});
    daily_weather = outerjoin(daily_weather, cloud_daily, 'Keys', 'date', 'MergeKeys', true);
end

% --- Wind speed ---
if ~isempty(wind_hourly)
    wind_hourly.date = dateshift(wind_hourly.time, 'start', 'day');
    [ag, ud] = findgroups(wind_hourly.date);
    wind_daily = table(ud, ...
        splitapply(@mean, wind_hourly.wind_speed_m_s, ag), ...
        'VariableNames', {'date', 'wind_speed_mean'});
    daily_weather = outerjoin(daily_weather, wind_daily, 'Keys', 'date', 'MergeKeys', true);
end

fprintf('  Weather daily table: %d rows, %d columns\n', height(daily_weather), width(daily_weather));

%% ============================================================
%% 7. MERGE AND COMPUTE CORRELATIONS
%% ============================================================

fprintf('\nMerging and computing correlations...\n');

% Ensure consistent datetime format (no timezone)
daily_network.date.TimeZone = '';
daily_weather.date.TimeZone = '';

merged = innerjoin(daily_network, daily_weather, 'Keys', 'date');
fprintf('  Merged table: %d rows\n', height(merged));

% Remove rows with NaN in critical columns
critical_cols = {'mean_power_kW', 'T_air_mean', 'T_air_max'};
valid_rows = true(height(merged), 1);
for c = 1:length(critical_cols)
    if ismember(critical_cols{c}, merged.Properties.VariableNames)
        valid_rows = valid_rows & isfinite(merged.(critical_cols{c}));
    end
end
merged = merged(valid_rows, :);
fprintf('  After NaN removal: %d rows\n', height(merged));

% --- Build list of available variables for correlation ---
demand_vars = {'total_energy_kWh', 'mean_power_kW', 'mean_std_power_kW', 'cv_power', ...
               'mean_active_hours', 'mean_flow_m3_h'};
weather_vars = {};
candidate_weather = {'T_air_mean', 'T_air_max', 'T_air_min', 'T_air_range', ...
                     'T_soil_mean', 'radiation_mean', 'radiation_max', 'radiation_total', ...
                     'cloud_cover_mean', 'wind_speed_mean'};
for v = 1:length(candidate_weather)
    if ismember(candidate_weather{v}, merged.Properties.VariableNames)
        weather_vars{end+1} = candidate_weather{v}; %#ok<SAGROW>
    end
end

all_vars = [demand_vars, weather_vars];
existing_vars = intersect(all_vars, merged.Properties.VariableNames, 'stable');

% Compute correlation matrix
corr_data_matrix = zeros(height(merged), length(existing_vars));
for v = 1:length(existing_vars)
    col = merged.(existing_vars{v});
    if ~isnumeric(col)
        col = str2double(string(col));
    end
    corr_data_matrix(:, v) = col;
end

corr_matrix = corr(corr_data_matrix, 'rows', 'pairwise');

fprintf('\n=== CORRELATION MATRIX ===\n');
fprintf('%25s', '');
for j = 1:length(existing_vars)
    fprintf('%15s', existing_vars{j}(1:min(14,end)));
end
fprintf('\n');
for i = 1:length(existing_vars)
    fprintf('%25s', existing_vars{i});
    for j = 1:length(existing_vars)
        fprintf('%15.3f', corr_matrix(i,j));
    end
    fprintf('\n');
end

%% ============================================================
%% 8. VISUALIZATION
%% ============================================================

fprintf('\nGenerating plots...\n');

output_folder = 'results/correlation_analysis';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% --- Figure 1: Demand vs Temperature scatter plots ---
fig1 = figure('Position', [100 100 1400 900], 'Name', 'Demand vs Temperature');

subplot(2,3,1);
scatter(merged.T_air_mean, merged.mean_power_kW, 8, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Daily Mean Air Temp [°C]'); ylabel('Mean Power [kW]');
r = corr(merged.T_air_mean, merged.mean_power_kW, 'rows', 'pairwise');
title(sprintf('Power vs T_{air,mean}\nr = %.3f', r));
grid on;

subplot(2,3,2);
scatter(merged.T_air_max, merged.mean_power_kW, 8, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Daily Max Air Temp [°C]'); ylabel('Mean Power [kW]');
r = corr(merged.T_air_max, merged.mean_power_kW, 'rows', 'pairwise');
title(sprintf('Power vs T_{air,max}\nr = %.3f', r));
grid on;

subplot(2,3,3);
if ismember('T_soil_mean', merged.Properties.VariableNames)
    scatter(merged.T_soil_mean, merged.mean_power_kW, 8, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel('Daily Mean Soil Temp [°C]'); ylabel('Mean Power [kW]');
    r = corr(merged.T_soil_mean, merged.mean_power_kW, 'rows', 'pairwise');
    title(sprintf('Power vs T_{soil}\nr = %.3f', r));
    grid on;
end

subplot(2,3,4);
scatter(merged.T_air_mean, merged.cv_power, 8, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Daily Mean Air Temp [°C]'); ylabel('CV of Power [-]');
r = corr(merged.T_air_mean, merged.cv_power, 'rows', 'pairwise');
title(sprintf('Power Stability vs T_{air,mean}\nr = %.3f', r));
grid on;

subplot(2,3,5);
scatter(merged.T_air_max, merged.cv_power, 8, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Daily Max Air Temp [°C]'); ylabel('CV of Power [-]');
r = corr(merged.T_air_max, merged.cv_power, 'rows', 'pairwise');
title(sprintf('Power Stability vs T_{air,max}\nr = %.3f', r));
grid on;

subplot(2,3,6);
if ismember('radiation_total', merged.Properties.VariableNames)
    scatter(merged.radiation_total, merged.mean_power_kW, 8, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel('Daily Total Radiation [Wh/m²]'); ylabel('Mean Power [kW]');
    r = corr(merged.radiation_total, merged.mean_power_kW, 'rows', 'pairwise');
    title(sprintf('Power vs Radiation\nr = %.3f', r));
    grid on;
end

saveas(fig1, fullfile(output_folder, 'demand_vs_temperature.png'));

% --- Figure 2: Time series overview ---
fig2 = figure('Position', [100 100 1400 1000], 'Name', 'Time Series Overview');

subplot(4,1,1);
plot(merged.date, merged.mean_power_kW, 'b-', 'LineWidth', 0.5);
ylabel('Mean Power [kW]'); title('Network Mean Heating Power');
grid on; xlim([t_start t_end]);

subplot(4,1,2);
plot(merged.date, merged.T_air_mean, 'r-', 'LineWidth', 0.5); hold on;
plot(merged.date, merged.T_air_max, 'r--', 'LineWidth', 0.5);
legend_entries = {'T_{air,mean}', 'T_{air,max}'};
if ismember('T_soil_mean', merged.Properties.VariableNames)
    plot(merged.date, merged.T_soil_mean, 'b-', 'LineWidth', 0.5);
    legend_entries{end+1} = 'T_{soil}';
end
legend(legend_entries, 'Location', 'best');
ylabel('[°C]'); title('Temperatures');
grid on; xlim([t_start t_end]);

subplot(4,1,3);
plot(merged.date, merged.cv_power, 'k-', 'LineWidth', 0.5);
ylabel('CV [-]'); title('Coefficient of Variation of Hourly Power (stability indicator)');
grid on; xlim([t_start t_end]);

subplot(4,1,4);
plot(merged.date, merged.mean_flow_m3_h, 'g-', 'LineWidth', 0.5);
ylabel('Flow [m³/h]'); title('Network Mean Flow Rate');
grid on; xlim([t_start t_end]);

saveas(fig2, fullfile(output_folder, 'time_series_overview.png'));

% --- Figure 3: Threshold analysis for T_air_max ---
fig3 = figure('Position', [100 100 1200 800], 'Name', 'Threshold Analysis');

thresholds = 0:1:25;
mean_power_at_threshold = nan(size(thresholds));
cv_at_threshold = nan(size(thresholds));
active_hours_at_threshold = nan(size(thresholds));

for ti = 1:length(thresholds)
    mask = merged.T_air_max < thresholds(ti);
    if sum(mask) > 10
        mean_power_at_threshold(ti) = mean(merged.mean_power_kW(mask));
        cv_at_threshold(ti) = mean(merged.cv_power(mask));
        active_hours_at_threshold(ti) = mean(merged.mean_active_hours(mask));
    end
end

subplot(2,2,1);
yyaxis left
plot(thresholds, mean_power_at_threshold, 'b-o', 'LineWidth', 1.5);
ylabel('Mean Power [kW]');
yyaxis right
n_days = arrayfun(@(th) sum(merged.T_air_max < th), thresholds);
plot(thresholds, n_days, 'r--', 'LineWidth', 1.5);
ylabel('Number of Days');
xlabel('T_{air,max} threshold [°C]');
title('Mean Power & Available Days vs Threshold');
grid on; legend('Mean Power', 'N Days', 'Location', 'best');

subplot(2,2,2);
plot(thresholds, cv_at_threshold, 'k-o', 'LineWidth', 1.5);
xlabel('T_{air,max} threshold [°C]');
ylabel('Mean CV of Power [-]');
title('Power Stability vs T_{air,max} Threshold');
grid on;
xline(10, 'r--', 'LineWidth', 2, 'Label', 'Proposed: 10°C');

subplot(2,2,3);
plot(thresholds, active_hours_at_threshold, 'g-o', 'LineWidth', 1.5);
xlabel('T_{air,max} threshold [°C]');
ylabel('Mean Active Hours per House');
title('Active Hours vs T_{air,max} Threshold');
grid on;
xline(10, 'r--', 'LineWidth', 2, 'Label', 'Proposed: 10°C');

subplot(2,2,4);
plot(merged.date, merged.T_air_max, 'b.', 'MarkerSize', 3); hold on;
yline(10, 'r-', 'LineWidth', 2, 'Label', 'T_{max} = 10°C');
yline(12, 'g--', 'LineWidth', 1, 'Label', 'T_{max} = 12°C');
yline(8, 'm--', 'LineWidth', 1, 'Label', 'T_{max} = 8°C');
ylabel('Daily Max Air Temp [°C]');
title('Daily T_{air,max} with Candidate Thresholds');
grid on; xlim([t_start t_end]);

saveas(fig3, fullfile(output_folder, 'threshold_analysis.png'));

% --- Figure 4: Radiation impact ---
if ismember('radiation_total', merged.Properties.VariableNames)
    fig4 = figure('Position', [100 100 1200 500], 'Name', 'Radiation Analysis');
    
    subplot(1,2,1);
    month_num = month(merged.date);
    scatter(merged.radiation_total, merged.cv_power, 12, month_num, 'filled', 'MarkerFaceAlpha', 0.5);
    cb = colorbar; ylabel(cb, 'Month');
    colormap(hsv(12));
    caxis([1 12]);
    xlabel('Daily Total Radiation [Wh/m²]'); ylabel('CV of Power [-]');
    title('Power Stability vs Radiation (colored by month)');
    grid on;
    
    subplot(1,2,2);
    shoulder = merged(merged.T_air_max >= 8 & merged.T_air_max <= 14, :);
    if height(shoulder) > 20
        scatter(shoulder.radiation_total, shoulder.cv_power, 15, shoulder.T_air_max, 'filled');
        cb = colorbar; ylabel(cb, 'T_{air,max} [°C]');
        xlabel('Daily Total Radiation [Wh/m²]'); ylabel('CV of Power [-]');
        title(sprintf('Shoulder Season (T_{max} 8-14°C)\nn=%d days', height(shoulder)));
        grid on;
    end
    
    saveas(fig4, fullfile(output_folder, 'radiation_analysis.png'));
end

% --- Figure 5: Correlation heatmap ---
fig5 = figure('Position', [100 100 900 700], 'Name', 'Correlation Heatmap');
imagesc(corr_matrix);
colorbar;
colormap(parula);  % Use parula (available in all MATLAB versions)
caxis([-1 1]);
xticks(1:length(existing_vars));
yticks(1:length(existing_vars));
short_names = strrep(existing_vars, '_', ' ');
xticklabels(short_names);
yticklabels(short_names);
xtickangle(45);
title('Correlation Matrix: Demand & Weather Variables');

% Add correlation values as text
for i = 1:length(existing_vars)
    for j = 1:length(existing_vars)
        text_color = 'k';
        if abs(corr_matrix(i,j)) > 0.7
            text_color = 'w';
        end
        text(j, i, sprintf('%.2f', corr_matrix(i,j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', text_color);
    end
end

saveas(fig5, fullfile(output_folder, 'correlation_heatmap.png'));

%% ============================================================
%% 9. SUMMARY STATISTICS
%% ============================================================

fprintf('\n=== KEY FINDINGS ===\n\n');

fprintf('--- Correlation with Mean Power (heating demand) ---\n');
for v = 1:length(weather_vars)
    vname = weather_vars{v};
    r = corr(merged.mean_power_kW, merged.(vname), 'rows', 'pairwise');
    fprintf('  %-20s : r = %+.3f\n', vname, r);
end

fprintf('\n--- Correlation with Power CV (instability indicator) ---\n');
for v = 1:length(weather_vars)
    vname = weather_vars{v};
    r = corr(merged.cv_power, merged.(vname), 'rows', 'pairwise');
    fprintf('  %-20s : r = %+.3f\n', vname, r);
end

fprintf('\n--- Threshold Comparison (T_air_max) ---\n');
fprintf('  %5s  %5s  %6s  %10s  %6s  %10s\n', ...
    'Tmax<', 'Days', '%', 'MeanPwr', 'CV', 'ActiveH');
for th = [6, 8, 10, 12, 14, 16, 18, 20]
    mask = merged.T_air_max < th;
    n = sum(mask);
    if n > 10
        fprintf('  %3d°C  %5d  %5.1f%%  %8.2f kW  %5.2f  %8.1f h\n', ...
            th, n, 100*n/height(merged), ...
            mean(merged.mean_power_kW(mask)), ...
            mean(merged.cv_power(mask)), ...
            mean(merged.mean_active_hours(mask)));
    end
end

fprintf('\n--- Seasonal Breakdown ---\n');
seasons = {'Winter (DJF)', 'Spring (MAM)', 'Summer (JJA)', 'Autumn (SON)'};
season_months = {[12,1,2], [3,4,5], [6,7,8], [9,10,11]};
fprintf('  %-15s  %10s  %6s  %10s  %6s\n', 'Season', 'MeanPwr', 'CV', 'T_air_max', 'Days');
for s = 1:4
    mask = ismember(month(merged.date), season_months{s});
    if sum(mask) > 10
        fprintf('  %-15s  %8.2f kW  %5.2f  %8.1f°C  %5d\n', ...
            seasons{s}, ...
            mean(merged.mean_power_kW(mask)), ...
            mean(merged.cv_power(mask)), ...
            mean(merged.T_air_max(mask)), ...
            sum(mask));
    end
end

% Save merged data for further analysis
writetable(merged, fullfile(output_folder, 'merged_daily_data.csv'));
fprintf('\nMerged daily data saved to %s\n', fullfile(output_folder, 'merged_daily_data.csv'));

fprintf('\n=== Analysis complete. Plots saved to %s ===\n', output_folder);



%% ============================================================
%% 2. HELPER FUNCTION: LOAD WEATHER FILE
%% ============================================================

 function weather_table = load_weather_file(filepath, value_col_name)
        % Loads a weather CSV with columns: date, value
        % Returns a table with columns: time (datetime), <value_col_name> (double)
        
        raw = readtable(filepath, 'TextType', 'string', 'Delimiter', ',');
        
        % The first column is the timestamp, second is the value
        time_strings = raw{:, 1};
        value_raw = raw{:, 2};
        
        % Convert values to double
        if isnumeric(value_raw)
            values = value_raw;
        elseif isstring(value_raw) || iscellstr(value_raw)
            values = str2double(value_raw);
        else
            values = double(value_raw);
        end
        
        % Strip timezone suffix before parsing
        % Handles formats like "2020-12-31 23:00:00+00:00"
        time_strings = regexprep(time_strings, '[+-]\d{2}:\d{2}$', '');
        time_strings = strtrim(time_strings);
        
        % Parse timestamps
        try
            times = datetime(time_strings, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        catch
            try
                times = datetime(time_strings, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
            catch
                % Last resort: let MATLAB figure it out
                fprintf('  WARNING: Trying automatic datetime parsing for %s\n', value_col_name);
                times = datetime(time_strings);
            end
        end
        
        % Ensure no timezone attached
        if ~isempty(times.TimeZone)
            times.TimeZone = '';
        end
        
        % Build output table
        weather_table = table(times, values, 'VariableNames', {'time', value_col_name});
        
        % Remove NaN values
        weather_table = weather_table(isfinite(weather_table.(value_col_name)), :);
        weather_table = sortrows(weather_table, 'time');
        
        fprintf('  %-25s: %d valid records (%.1f%% of raw)\n', ...
            value_col_name, height(weather_table), ...
            100 * height(weather_table) / length(values));
    end
