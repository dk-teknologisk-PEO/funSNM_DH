function [meter_data, network_temp, topology] = importData(config,network)

city = config.project.location.city;
data_path = config.project.paths.data;
topology_path = config.project.paths.topology;
start_time = config.project.time.start;
end_time = config.project.time.end;

if config.project.location.city == "Aalborg"
    [meter_data, network_temp, topology] = importAalborgData(data_path,topology_path, network, "hourly");
    meter_data = meter_data(meter_data.timestamp>=start_time,:);
    meter_data = meter_data(meter_data.timestamp<=end_time,:);
    
    meter_data.flow_kg_h = round(meter_data.flow_kg_s.*3600,5);
    meter_data.T_supply_C = round(meter_data.T_supply_C,2);
    meter_data.timestamp = split(string(meter_data.timestamp),'+');
    meter_data.timestamp = meter_data.timestamp(:,1);
    meter_data.timestamp = datetime(meter_data.timestamp);
    meter_data = meter_data(meter_data.timestamp.Minute==00,:);
    meter_data = meter_data((meter_data.timestamp.Hour>23) | (meter_data.timestamp.Hour<4),:);
    meter_data = meter_data(:,{'timestamp', 'house_id','T_supply_C', 'flow_kg_h','T_tapoff_C'});

    network_temp = network_temp(network_temp.timestamp>=start_time,:);
    network_temp = network_temp(network_temp.timestamp<=end_time,:);
end


