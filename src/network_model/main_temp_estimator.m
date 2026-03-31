function T_main_sq = main_temp_estimator(parameters, current_data, current_data_sub, current_T_soil_C, U_main, T_start_C)
% MAIN_TEMP_ESTIMATOR Cost function for fminsearch.
% Computes the sum of squared differences between the forward-simulated
% main pipe temperature profile and the back-calculated T_main from meters.
%
% parameters(1) = T_start_C (inlet temperature to optimize)    
    
    if nargin == 5 || isempty(T_start_C) || isnan(T_start_C)
        T_start_C = parameters;
    end
    
        
    % Segment lengths along main pipe
    xDiff = current_data.x_pos_m;
    xDiff(2:end) = diff(xDiff);
    
    % Cumulative flow (total flow in each main pipe segment)
    Q_main_total = cumsum(current_data.flow_kg_h, 1, 'reverse');
    
    % move expected temperatures in the main pipe (from utility meter
    % measurements) into seperate vector for better readability 
    temperatures = current_data_sub.T_supply_C;
    
    % Forward simulate temperature along main pipe
    num_houses = size(current_data, 1);
    temp = nan(num_houses, 1);
    dT = zeros(num_houses, 1);
    T_prev = T_start_C;

    for t = 1:num_houses
        temp(t) = get_supply_temp(T_prev, Q_main_total(t), U_main, xDiff(t), current_T_soil_C);
        house_id = current_data.house_id(t);

        if ismember(house_id, current_data_sub.house_id)
            ref_T_main = current_data.T_main_C(t);
            if isfinite(ref_T_main)
                dT(t) = temp(t) - ref_T_main;
            end
        end
        T_prev = temp(t);
    end

    T_main_sq = dT' * dT;
end