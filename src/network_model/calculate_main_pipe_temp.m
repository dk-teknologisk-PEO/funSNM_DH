function T_main_fit_C = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_cutoff, states, T_start_C)

    
    Q_total = cumsum(current_data.flow_kg_h,'reverse');
    num_houses = size(current_data,1);
    T_main_C = nan([num_houses,1]);
    x_diff_m = [current_data.x_pos_m];
    x_diff_m(2:end) = diff(x_diff_m);
    T_main_fit_C = nan([num_houses,1]);

    if nargin == 6
        T_prev_C = T_start_C;
    
        for i=1:size(current_data,1)
            T_prev_C = get_supply_temp(T_prev_C,Q_total(i),x_main(2),x_diff(i),T_ambient);
            T_main_C(i)=T_prev_C;
        end

    else
        T_start_C = NaN;
        u_service = cellfun(@(s) s.x(2), states);
        offsets = cellfun(@(s) s.x(1), states);
        service_lengths_m = [current_data.length_service_m];
        T_main_C = get_main_temp(current_data.T_supply_C+offsets', current_data.flow_kg_h,u_service', service_lengths_m, current_T_soil_C);
        current_data.T_main_C = T_main_C;
        T_first_guess = max(current_data.T_main_C);
        for m = 1:num_houses
        
            current_data_sub = current_data;
            if current_data_sub(m,:).flow_kg_h < flow_cutoff
                continue
            else
                current_id = current_data(m,:).house_id;
                current_data_sub(m,:)=[];
                fun = @(x)main_temp_estimator(x, current_data,current_data_sub, current_T_soil_C, U_csac, T_start_C);
                p = fminsearch(fun,T_first_guess);
                T_first_guess = p(1);
                T_prev_C = p(1);
                
                for i=1:m
                    T_prev_C = get_supply_temp(T_prev_C,Q_total(i),U_csac,x_diff_m(i),current_T_soil_C);
                    if i==m
                        T_main_fit_C(m) = T_prev_C;
                    end
                end
            end
        end
    end
