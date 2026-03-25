function [T_soil_C, T_air_C] = soilTemp(config)

api_key = '912b1b04-ff7a-4383-b8c3-3fc0d9fa30e5';
city = config.project.location.city;
t_start = datetime(config.project.time.start);
t_end = datetime(config.project.time.end);
start_path = config.project.paths.weather;
path_soil = strcat(start_path,city,"_soil.mat");
path_air = strcat(start_path, city, "_air.mat");

try
    load(path_soil,"T_soil_C");
catch
    T_soil_C = get_soil_temperature(city,t_start,t_end,api_key);
    save(path_soil,"T_soil_C")
end


try
    load(path_air,"T_air_C");
catch
    T_air_C = get_air_temperature(city,t_start,t_end,api_key);
    save(path_air,"T_air_C");
end