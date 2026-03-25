function T_air = get_air_temperature(city, t_start,t_end, api_key)%lat, lon, number_of_weather_stations, combined_data)
% function for finding the soil temperature (at 30 cm below ground) for a
% given city. The function takes city name (might have issues if multiple
% cities have the same name), finds the longitude and lattitude, finds the
% 10km-grid ID, and then uses DMI to find the soil temperature in the
% period from t_start to t_end.

    
    % find longitude and latitude of the city
    [long, lat] = city_position(city);
    % get 10km-grid id
    grid_10km = get_10km_grid(long,lat);
    % api_key = '912b1b04-ff7a-4383-b8c3-3fc0d9fa30e5';
        
    % Generate all months between start and end dates
    dates = t_start:calmonths(1):t_end;
    
    for date=dates
        date_interval = strcat(string(year(date)), '-', sprintf('%02d',month(date)), '-01T00:00:00Z/',string(year(date+calmonths(1))),'-', sprintf('%02d',month(date+calmonths(1))),'-01T00:00:00Z');
        disp(date_interval)

        % use the DMI API to find the data
        data = webread('https://dmigw.govcloud.dk/v2/climateData/collections/10kmGridValue/items?','cellId',grid_10km,'parameterId','mean_temp','datetime',date_interval,'api-key',api_key);
        % data = webread('https://dmigw.govcloud.dk/v2/climateData/collections/10kmGridValue/items?','cellId',grid_10km,'parameterId','temp_soil_30','datetime',date_interval,'api-key',api_key);
        % initialize vectors to contain soil temperature and timestamps
        values = zeros(size(data.features,1),1);
        time = NaT(size(data.features,1),1);

        % run through all data to do sanitation. Only use hourly values,
        % not averages over day/month/year
        for j=1:length(values)
            if ~strcmp(data.features(j).properties.timeResolution,'day')
                continue
            end
            % extract timestamp. Use end of hour, as this corresponds to
            % the time in which we have data from the district heating
            % meters
            timestamp = data.features(j).properties.to;
            time_string = strsplit(timestamp,'+');
            % check if it is summer time, and correct if required
            offset_time = strsplit(time_string{2},':');
            offset_hour = str2double(offset_time{1});
            utc_time = datetime(time_string{1}) + hours(offset_hour);
            time(j) = utc_time;
            % extract the soil temperature
            values(j) = data.features(j).properties.value;
        end
        T_air_temporary = rmmissing(timetable(time,values));

        % add to main timetable
        if ~exist("T_air")
            T_air = T_air_temporary;
        else
            T_air = [T_air;T_air_temporary;];
        end
    end
    % remove rows with no timestamp
    T_air(isnat(T_air.time),:)=[];
    % sort the table chronologically
    T_air = sortrows(T_air,"time", "ascend");
    % remove duplicates
    T_air = unique(T_air);
    % remove bad data (soil temperature shouldn't change more than 2 °C/h
    % to avoid issues with first value being bad, use the average of the
    % first day as first value
    % T_soil.values(1) = mean(T_soil.values(1:24));
    % for t = 2:height(T_soil)
    %     if abs(T_soil.values(t)-T_soil.values(t-1))>2
    %         T_soil.values(t) = T_soil.values(t-1);
    %     end
    % end
    % take a rolling mean.
    % T_soil.values = movmean(T_soil.values,10);
end