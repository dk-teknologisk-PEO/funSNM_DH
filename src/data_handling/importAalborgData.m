function [house_meter_data_long, network_temp, topology] = importAalborgData(path_data, path_topology, network_ID, freq)

% path_data = "C:\Users\peo\Documents\GitHub\DH simulation\results";
% path_topology = "C:\Users\peo\Documents\GitHub\DH simulation\network_topologies";

w = waitbar(0, "please wait");

network_ID = sprintf('%02d',network_ID);

topology = fileread(strcat(path_topology, "\network_", network_ID, ".json"));
topology = jsondecode(topology);

data_file = strcat(path_data, "\network_", network_ID,"\", freq, "_consumer_data.csv");
network_file = strcat(path_data, "\network_", network_ID,"\", freq, "_network_temperatures.csv");

try
    house_meter_data_long = readtable(data_file);
    network_temp = readtable(network_file);
catch
    
    temp_file = strcat(path_data, "\network_", network_ID,"\", freq, "_temperatures.csv");
    flow_file = strcat(path_data, "\network_", network_ID,"\", freq, "_flows.csv");
    
    temp_data = readtable(temp_file, "Delimiter", ',');
    temp_data.timestamp = datetime(temp_data.timestamp,"InputFormat","uuuu-MM-dd HH:mm:ss+00:00");
    waitbar(1/100, w)
    
    flow_data = readtable(flow_file, "Delimiter", ',');
    flow_data.timestamp = datetime(flow_data.timestamp,"InputFormat","uuuu-MM-dd HH:mm:ss+00:00");
    waitbar(2/100, w)
    
    full_data = join(temp_data, flow_data);
    
    house_ids = [topology.houses.id];
    time = full_data.timestamp;
    
    for t = 1:length(time)
        flow = NaN([length(house_ids),1]);
        temp = NaN([length(house_ids),1]);
        tapoff = NaN([length(house_ids),1]);
        
        timestamp = NaT([length(house_ids),1]);
        timestamp(:) = time(t);
        
        data_temp = full_data(full_data.timestamp == time(t),:);
        for h = 1:length(house_ids)
            house_id = house_ids(h);
            csac = topology.houses(h).cul_de_sac_id;
            temp_col = strcat('H_',string(house_id),'_s');
            tapoff_col = strcat('J_CSAC',string(csac),'_H',string(house_id),'_s');
            flow_col = strcat('service_H', string(house_id));
            temp(h) = table2array(data_temp(:,temp_col));
            tapoff(h) = table2array(data_temp(:,tapoff_col));
            flow(h) = table2array(data_temp(:,flow_col));
        end
    
        table_temp = table(timestamp, house_ids',temp, tapoff,flow);
        table_temp.Properties.VariableNames = {'timestamp','house_id','T_supply_C','T_tapoff_C','flow_kg_s'};
        try
            house_meter_data_long = [house_meter_data_long;table_temp];
            if mod(t,100)==0
                waitbar((t+2)/length(time),w)
            end
        catch
            house_meter_data_long = table_temp;
        end
    end
    waitbar(3/100, w)
    
    network_temp = temp_data;
    for h = 1:length(house_ids)
        house_id = house_ids(h);
        csac = topology.houses(h).cul_de_sac_id;
        temp_col = strcat('H_',string(house_id),'_s');
        tapoff_col = strcat('J_CSAC',string(csac),'_H',string(house_id),'_s');
        network_temp(:,temp_col) = [];
        network_temp(:,tapoff_col) = [];
    end
    writetable(house_meter_data_long, data_file)
    writetable(network_temp, network_file)
end
close(w)