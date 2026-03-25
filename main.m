% Clear workspace, add paths, set parameters
clear all
close all
clc

addpath('src/kalman_filter', 'src/network_model', 'src/data_handling', 'src/diagnostics', 'src/particle_filter', 'config')

% read config-file
config = jsondecode(fileread("config.json"));

% ukf-configuration
R_base = config.project.initialization.ukf.measurement_noise^2;
Q_base = diag([(config.project.initialization.ukf.process_noise_offset)^2, (config.project.initialization.ukf.process_noise_U)^2]);
P_base = diag([(config.project.initialization.ukf.state_uncertainty_offset)^2, (config.project.initialization.ukf.state_uncertainty_U)^2]);

% pf-configuration
num_particles = config.project.initialization.pf.num_particles;

% gating configuration
flow_threshold = config.project.cutoff.flow_cutoff;
delta_T_gate_threshold = config.project.cutoff.delta_T_gate_threshold;
min_active_houses = config.project.initialization.min_active_houses;
max_innovation = config.project.initialization.max_innovation_C;
air_temp_cutoff = config.project.initialization.max_air_temperature;

networks = config.project.datasets.datasets;

% load network topology and measurement data
[T_soil_C, T_air_C] = soilTemp(config);
output_folder_ukf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'),'/ukf'); % Create a unique folder for this run
output_folder_pf = fullfile('results', datestr(now, 'yyyy-mm-dd_HHMM'),'/pf'); % Create a unique folder for this run
w = waitbar(0.0, "Starting analysis");


% run through all the networks from the config-file
for network = networks
    % load relevant data for the network
    [meter_data, network_data, topology] = importData(config,network);
    timestamps = unique(meter_data.timestamp);
    
    % find the ids of the csacs from the topology
    csac_ids = [topology.cul_de_sacs.id];
    prev_houses = 0;

    % loop through all cul-de-sacs in the network
    for csac=csac_ids
        
        % extract relevant data for the csac
        U_csac = topology.pipe_parameters.csac_pipe.insulation_W_m_K; % U-value of the main pipe
        houses_csac = ([topology.houses.cul_de_sac_id]==csac);
        topology_csac = topology.houses(houses_csac);
        num_houses_csac = sum(houses_csac);
        houses_csac_ids = [topology.houses(houses_csac).id];
        true_offset = zeros([size(topology_csac),1]);

        % create a table with reference-values (only possible when there is
        % a topology)
        ground_truth_csac = table([topology_csac.id]', true_offset, [topology_csac.service_pipe_insulation_W_m_K]', ...
            'VariableNames', {'house_id', 'true_offset', 'true_U'});
        
        % create a table with lengths and positions of service pipes.
        csac_table = table([topology_csac.id]', [topology_csac.service_pipe_len_m]', [topology_csac.dist_on_cul_de_sac_m]', 'VariableNames',{'house_id','length_service_m','x_pos_m'});


        % extract relevant meter data
        meter_data_csac = meter_data(ismember(meter_data.house_id,houses_csac_ids),:);
        meter_data_csac = join(meter_data_csac, csac_table, "Keys", "house_id");

        % initialize a Kalman filter and particle filter for each house on the csac and store
        % them i a structure for easy access
        ukf_states = cell(1, num_houses_csac);
        pf_particles = cell(1, num_houses_csac);
        pf_states = cell(1, num_houses_csac);
        for i = 1:num_houses_csac
            x_init = [randn()*0.3; 0.14 + randn()*0.03]; % make initial guess for house i
            
            % set up the UKF-state for house i
            ukf_states{i}=initialize_kalman_filter(x_init, Q_base, R_base, P_base);
            % set up particles for house i
            pf_particles{i} = initialize_pf_state(num_particles, x_init, P_base);
            pf_states{i}.x = x_init;
            pf_states{i}.P = P_base;
        end
        
        % set up loggers for the two filters
        logger_ukf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);
        logger_pf = initialize_logger(num_houses_csac, length(timestamps), houses_csac_ids);

        % keep track of the overall offsets, allowing us to modify their
        % means later on
        ukf_offsets = nan([num_houses_csac,1]);
        pf_offsets = nan([num_houses_csac,1]);

        % run through the timesteps
        for t = 1:length(timestamps)
            disp(t)
            if mod(t,100)==0
                waitbar(t/(length(timestamps)*length(csac_ids))+(csac)/length(csac_ids),w,strcat("Running network. At timestep: ", string(t), "/", string(length(timestamps)),", csac: ",string(csac),"/",string(length(csac_ids)),", network: ", string(network)))
            end
            time = timestamps(t);

            % extract data from this timestep
            current_data = meter_data_csac(meter_data_csac.timestamp==time,:);
            current_T_soil_C = T_soil_C(T_soil_C.time==time,:).values;
            current_T_air_C = T_air_C((year(T_air_C.time)==year(time)) & (month(T_air_C.time)==month(time)) & (day(T_air_C.time)==day(time)),:).values;
            num_active_houses = sum(current_data.flow_kg_h >= flow_threshold);
            is_csac_active = (num_active_houses > config.project.initialization.min_active_houses);
            skip_timestep = false;
            if isempty(current_data) || isempty(current_T_soil_C) || (current_T_air_C > air_temp_cutoff )
                skip_timestep = true;
            else
                T_junction_ukf_C = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_threshold, ukf_states);
                current_data.T_main_ukf_C = T_junction_ukf_C;
    
                T_junction_pf_C = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_threshold, pf_states);
                current_data.T_main_pf_C = T_junction_pf_C;
            end

            for i = 1:num_houses_csac
                log_ukf = false;
                log_pf = false;
                if ~skip_timestep
                    house_id = houses_csac_ids(i);
                    house_data = current_data(current_data.house_id==house_id,:);
                    is_flow_sufficient = (house_data.flow_kg_h > config.project.cutoff.flow_cutoff);
                   

    
                    delta_T_ukf_sufficient = ((house_data.T_main_ukf_C - current_T_soil_C) >= delta_T_gate_threshold);
                    delta_T_pf_sufficient = ((house_data.T_main_pf_C - current_T_soil_C) >= delta_T_gate_threshold);
                    if delta_T_ukf_sufficient && is_flow_sufficient && is_csac_active
                    % update UKF and particle filter
                        predicted_temp = get_supply_temp(house_data.T_main_ukf_C, house_data.flow_kg_h,ukf_states{i}.x(2), house_data.length_service_m, current_T_soil_C);
                        innovation = house_data.T_supply_C - (predicted_temp - ukf_states{i}.x(1));
                        if abs(innovation) < max_innovation
                            [ukf_states{i}, diagnostics_ukf] = update_ukf_house(ukf_states{i}, house_data, current_T_soil_C);
                            ukf_offsets(i) = ukf_states{i}.x(1);
                            log_ukf = true;
                        end
                    end
                    
                    if delta_T_pf_sufficient && is_flow_sufficient && is_csac_active
                        predicted_temp = get_supply_temp(house_data.T_main_pf_C, house_data.flow_kg_h,pf_states{i}.x(2), house_data.length_service_m, current_T_soil_C);
                        innovation = house_data.T_supply_C - (predicted_temp - pf_states{i}.x(1));
                        if abs(innovation) < max_innovation
                            [pf_particles{i}, est_pf, cov_pf, diagnostics_pf] = update_pf_house(pf_particles{i}, house_data, current_T_soil_C, R_base, Q_base);
                            pf_states{i}.x = est_pf;
                            pf_states{i}.P = cov_pf;
                            pf_offsets(i) = pf_states{i}.x(1);
                            log_pf = true;
                        end
                    end
                end
                if log_ukf
                    logger_ukf = update_logger(logger_ukf, t, i, time, ukf_states{i}, diagnostics_ukf);
                elseif t > 1
                    logger_ukf.state_estimates(:, i, t) = logger_ukf.state_estimates(:, i, t-1);
                    logger_ukf.covariance_posterior(:, i, t) = logger_ukf.covariance_posterior(:, i, t-1);
                end
                if log_pf
                    logger_pf = update_logger(logger_pf, t, i, time, pf_states{i}, diagnostics_pf);
                elseif t > 1
                    logger_pf.state_estimates(:, i, t) = logger_pf.state_estimates(:, i, t-1);
                    logger_pf.covariance_posterior(:, i, t) = logger_pf.covariance_posterior(:, i, t-1);

                end
            end
            mean_ukf_offsets = mean(ukf_offsets, 'omitmissing');
            mean_pf_offsets = mean(pf_offsets, 'omitmissing');
            for i=1:num_houses_csac
                if (~isnan(mean_ukf_offsets) && ~isnan(ukf_states{i}.x(1)))
                    ukf_states{i}.x(1)=ukf_states{i}.x(1)-mean_ukf_offsets;
                end
                if (~isnan(mean_pf_offsets) && ~isnan(pf_states{i}.x(1)))
                    pf_states{i}.x(1)=pf_states{i}.x(1)-mean_pf_offsets;
                    pf_particles{i}(1,:) = pf_particles{i}(1,:) - mean_pf_offsets;
                end
            end
        end
        % prev_houses = prev_houses + size(num_houses_csac);
        % --- Plot and Save Diagnostics for this CSAC ---
        plot_diagnostics(logger_ukf, ground_truth_csac, csac, network, output_folder_ukf);
        plot_diagnostics(logger_pf, ground_truth_csac, csac, network, output_folder_pf)
    end % End of CSAC loop
end % End of network loop

close(w);
disp('Analysis complete. All diagnostic plots have been saved.')
