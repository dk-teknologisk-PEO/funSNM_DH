% function for finding the "goodness" of the temperature profile along the
% main pipe

function T_main_sq = main_temp_estimator(parameters, current_data,current_data_sub, current_T_soil_C, U_main, T_start_C)
    
    
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
    for t = 1:length(temp)
        house_id = current_data(t,:).house_id;
        temp(t) = get_supply_temp(T_prev,Q_main_total(t),U_main,xDiff(t),current_T_soil_C);
        if ismember(house_id,current_data_sub.house_id) %|| ((t==length(temp)) & (~isnan(T_end_measured))) %   
            dT(t) = temp(t)-current_data(t,:).T_main_C; 
            % temperatures(1)=[];
        end
        T_prev = temp(t);
    end
    T_main_sq = dT'*dT;
