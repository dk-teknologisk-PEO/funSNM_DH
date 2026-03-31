% new file: main_temp_estimator_global.m
function T_main_sq = main_temp_estimator_global(parameters, current_data, ref_indices, current_T_soil_C, U_main)
% MAIN_TEMP_ESTIMATOR Cost function for finding global T_start and master_offset.
    
    T_start_C = parameters(1);
    master_offset = parameters(2);
        
    xDiff = current_data.x_pos_m;
    xDiff(2:end) = diff(xDiff);
    Q_main_total = cumsum(current_data.flow_kg_h, 1, 'reverse');
    
    num_houses = size(current_data, 1);
    dT_vector = []; % Dynamically build vector of errors
    
    T_prev = T_start_C;
    for t = 1:num_houses
        temp = get_supply_temp(T_prev, Q_main_total(t), U_main, xDiff(t), current_T_soil_C);
        
        if ismember(t, ref_indices)
            ref_T_main = current_data.T_main_C(t) - master_offset;
            dT_vector(end+1) = temp - ref_T_main;
        end
        T_prev = temp;
    end

    T_main_sq = dT_vector * dT_vector'; % Sum of squares
end