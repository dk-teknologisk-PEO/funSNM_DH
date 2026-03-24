function [meter_data, network_temp, topology] = importData(config,network)

city = config.project.location.city;
data_path = config.project.paths.data;
topology_path = config.project.paths.topology;

if config.project.location.city == "Aalborg"
    [meter_data, network_temp, topology] = importAalborgData(data_path,topology_path, network, "hourly");

    meter_data.flow_kg_h = round(meter_data.flow_kg_s.*3600,5);
    meter_data.T_supply_C = round(meter_data.T_supply_C,2);
    meter_data.timestamp = split(string(meter_data.timestamp),'+');
    meter_data.timestamp = meter_data.timestamp(:,1);
    meter_data.timestamp = datetime(meter_data.timestamp);
    meter_data = meter_data(meter_data.timestamp.Minute==00,:);
    meter_data = meter_data((meter_data.timestamp.Hour>23) | (meter_data.timestamp.Hour<4),:);
    meter_data = meter_data(:,{'timestamp', 'house_id','T_supply_C', 'flow_kg_h','T_tapoff_C'});
end


