% function for finding the "goodness" of the temperature profile along the
% main pipe

function T_main_sq = main_temp_estimator(parameters, current_data,current_data_sub, current_T_soil_C, U_main, T_start_C)
    
<<<<<<< HEAD
    if nargin == 5 || isempty(T_start_C) || isnan(T_start_C)
        T_start_C = parameters(1);
        master_offset = parameters(2);
    end
=======
>>>>>>> parent of d5cf9ed (Improve temp models and filter stability)
    
    T_start_C = parameters;
    
    % move expected temperatures in the main pipe (from utility meter
    % measurements) into seperate vector for better readability 
    temperatures = current_data_sub.T_supply_C;
        
    % find distance between service pipe connections at main pipe
    xDiff = current_data.x_pos_m;
    xDiff(2:end) = diff(xDiff);
    
    % xDiff = nan(size(xPos));
    % xDiff(1) = xPos(1);
    % xDiff(2:end) = xPos(2:end)-xPos(1:end-1);

    
    if ~isnan(T_start_C)
        T_start_C = T_start_C;
    end

    % find flow in different segments of the main pipe
    Q_main_total=cumsum(current_data.flow_kg_h,1,"reverse");
        
    % find temperature at entrance to each subpipe
    temp = nan(size(current_data.flow_kg_h));
    dT = zeros(size(current_data.flow_kg_h));
    T_prev = T_start_C;
<<<<<<< HEAD

    for t = 1:num_houses
        temp(t) = get_supply_temp(T_prev, Q_main_total(t), U_main, xDiff(t), current_T_soil_C);
        house_id = current_data.house_id(t);

        if ismember(house_id, current_data_sub.house_id)
            % The reference T_main is back-calculated from the measured T_supply.
            % If all meters have a positive offset, the back-calculated T_main will be too high.
            % We adjust it by the master_offset we are trying to find.
            ref_T_main = current_data.T_main_C(t) - master_offset;
            if isfinite(ref_T_main)
                dT(t) = temp(t) - ref_T_main;
            end
=======
    for t = 1:length(temp)
        house_id = current_data(t,:).house_id;
        temp(t) = get_supply_temp(T_prev,Q_main_total(t),U_main,xDiff(t),current_T_soil_C);
        if ismember(house_id,current_data_sub.house_id) %|| ((t==length(temp)) & (~isnan(T_end_measured))) %   
            dT(t) = temp(t)-current_data(t,:).T_main_C; 
            % temperatures(1)=[];
>>>>>>> parent of d5cf9ed (Improve temp models and filter stability)
        end
        T_prev = temp(t);
    end
    T_main_sq = dT'*dT;
