% inspect_csac_summary.m
% Generates a compact per-CSAC summary of physical input conditions

clear all; close all; clc;
addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'config')

config = jsondecode(fileread("config.json"));
[meter_data, network_data, topology] = importData(config, 1);
[T_soil_C, T_air_C] = soilTemp(config);
timestamps = unique(meter_data.timestamp);

flow_cutoff = config.project.cutoff.flow_cutoff;
delta_T_threshold = config.project.cutoff.delta_T_gate_threshold;

csac_ids = [topology.cul_de_sacs.id];

for csac = csac_ids
    houses_idx = ([topology.houses.cul_de_sac_id] == csac);
    topology_csac = topology.houses(houses_idx);
    houses_ids = sort([topology_csac.id]);
    
    fprintf('\n========== CSAC %d ==========\n', csac);
    fprintf('Houses: %s\n', mat2str(houses_ids));
    
    % Print pipe lengths, positions, true U-values
    fprintf('\n  House | Pipe_L [m] | Dist_on_csac [m] | True_U [W/m/K]\n');
    fprintf('  ------|------------|------------------|---------------\n');
    for i = 1:length(topology_csac)
        h = topology_csac(i);
        fprintf('  %5d | %10.1f | %16.1f | %13.4f\n', ...
            h.id, h.service_pipe_len_m, h.dist_on_cul_de_sac_m, h.service_pipe_insulation_W_m_K);
    end
    
    % Per-house flow and delta-T statistics
    fprintf('\n  House | mean_flow | median_flow | pct_above_cutoff | mean_DeltaT | pct_DeltaT_ok | n_measurements\n');
    fprintf('  ------|-----------|-------------|------------------|-------------|---------------|---------------\n');
    
    for i = 1:length(houses_ids)
        hid = houses_ids(i);
        hdata = meter_data(meter_data.house_id == hid, :);
        
        % Match soil temp
        T_soil_matched = nan(height(hdata), 1);
        for j = 1:height(hdata)
            idx = find(T_soil_C.time == hdata.timestamp(j), 1);
            if ~isempty(idx)
                T_soil_matched(j) = T_soil_C.values(idx);
            end
        end
        
        % Estimate T_main using true U
        U_true = topology_csac(i).service_pipe_insulation_W_m_K;
        L_pipe = topology_csac(i).service_pipe_len_m;
        T_main_est = get_main_temp(hdata.T_supply_C, hdata.flow_kg_h, U_true, L_pipe, T_soil_matched);
        delta_T = T_main_est - T_soil_matched;
        
        flow = hdata.flow_kg_h;
        n = height(hdata);
        
        pct_flow_ok = 100 * sum(flow >= flow_cutoff) / n;
        pct_deltaT_ok = 100 * sum(delta_T >= delta_T_threshold, 'omitnan') / n;
        
        fprintf('  %5d | %9.1f | %11.1f | %15.1f%% | %11.1f | %12.1f%% | %14d\n', ...
            hid, mean(flow, 'omitnan'), median(flow, 'omitnan'), pct_flow_ok, ...
            mean(delta_T, 'omitnan'), pct_deltaT_ok, n);
    end
    
    % Monthly breakdown of flow activity for the whole CSAC
    fprintf('\n  Monthly flow activity (pct timesteps with >= %d houses above flow cutoff):\n', config.project.initialization.min_active_houses);
    months_to_check = 1:5;
    for m = months_to_check
        month_mask = month(meter_data.timestamp) == m & ismember(meter_data.house_id, houses_ids);
        month_data = meter_data(month_mask, :);
        month_ts = unique(month_data.timestamp);
        n_active_ts = 0;
        for t = 1:length(month_ts)
            ts_data = month_data(month_data.timestamp == month_ts(t), :);
            if sum(ts_data.flow_kg_h >= flow_cutoff) >= config.project.initialization.min_active_houses
                n_active_ts = n_active_ts + 1;
            end
        end
        fprintf('    Month %d: %d/%d timesteps active (%.1f%%)\n', ...
            m, n_active_ts, length(month_ts), 100*n_active_ts/max(1,length(month_ts)));
    end
end