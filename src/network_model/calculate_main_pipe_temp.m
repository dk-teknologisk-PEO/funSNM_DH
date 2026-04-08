
function [T_main_fit_C, master_offset_C] = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_cutoff, states, T_start_C)
% CALCULATE_MAIN_PIPE_TEMP Estimates main pipe temperature and a master offset.
% A single optimization is performed to find the best inlet temperature and
% a shared meter offset for the entire cul-de-sac.

    num_houses = size(current_data,1);
    
    % Initialize outputs
    T_main_fit_C = nan(num_houses,1);
    master_offset_C = NaN;
    
    % --- (This initial part of back-calculating T_main_C is still needed) ---
    u_service = cellfun(@(s) s.x(2), states(:));
    offsets   = cellfun(@(s) s.x(1), states(:));
    service_lengths_m = current_data.length_service_m;
    T_main_C = get_main_temp(current_data.T_supply_C + offsets, current_data.flow_kg_h, u_service, service_lengths_m, current_T_soil_C);
    T_main_C(~isfinite(T_main_C)) = NaN;
    current_data.T_main_C = T_main_C;
    
    % --- Global Optimization Step ---
    valid_houses_for_fit = find(isfinite(current_data.T_main_C));
    if numel(valid_houses_for_fit) < 2
        return; % Not enough data to do a fit
    end
    
    opts = optimset('Display', 'off', 'TolX', 0.05, 'TolFun', 0.05, 'MaxIter', 200);
    
    % Modify main_temp_estimator to accept a subset of indices to use
    fun = @(params) main_temp_estimator_global(params, current_data, valid_houses_for_fit, current_T_soil_C, U_csac);

    default_start_T = max(T_main_C(valid_houses_for_fit));
    start_vector = [default_start_T, 0];
    
    [p_global, ~, exitflag] = fminsearch(fun, start_vector, opts);
    
    if exitflag <= 0 || any(~isfinite(p_global))
        warning('Global T_main optimization failed.');
        return;
    end
    
    optimized_T_start = p_global(1);
    master_offset_C   = p_global(2);
    
    % --- Forward-simulate for ALL houses from the single optimized T_start ---
    Q_total = cumsum(current_data.flow_kg_h,'reverse');
    x_diff_m = [current_data.x_pos_m];
    x_diff_m(2:end) = diff(x_diff_m);

    T_prev_C = optimized_T_start;
    for i = 1:num_houses
        T_prev_C = get_supply_temp(T_prev_C, Q_total(i), U_csac, x_diff_m(i), current_T_soil_C);
        % We only assign a valid temperature if the house has flow
        if current_data.flow_kg_h(i) >= flow_cutoff
            T_main_fit_C(i) = T_prev_C;
        end
    end
end

% function T_main_fit_C = calculate_main_pipe_temp(current_data, current_T_soil_C, U_csac, flow_cutoff, states, T_start_C)
% 
% 
%     Q_total = cumsum(current_data.flow_kg_h,'reverse');
%     num_houses = size(current_data,1);
%     T_main_C = nan([num_houses,1]);
%     x_diff_m = [current_data.x_pos_m];
%     x_diff_m(2:end) = diff(x_diff_m);
%     T_main_fit_C = nan([num_houses,1]);
% 
%     if nargin == 6
%         T_prev_C = T_start_C;
% 
%         for i=1:size(current_data,1)
%             T_prev_C = get_supply_temp(T_prev_C,Q_total(i),x_main(2),x_diff(i),T_ambient);
%             T_main_C(i)=T_prev_C;
%         end
% 
%     else
%         T_start_C = NaN;
%         u_service = cellfun(@(s) s.x(2), states);
%         offsets = cellfun(@(s) s.x(1), states);
%         service_lengths_m = [current_data.length_service_m];
%         T_main_C = get_main_temp(current_data.T_supply_C+offsets', current_data.flow_kg_h,u_service', service_lengths_m, current_T_soil_C);
%         current_data.T_main_C = T_main_C;
% 
%         for m = 1:num_houses
% 
%             current_data_sub = current_data;
%             if current_data_sub(m,:).flow_kg_h < flow_cutoff
%                 continue
%             else
%                 current_id = current_data(m,:).house_id;
%                 current_data_sub(m,:)=[];
%                 fun = @(x)main_temp_estimator(x, current_data,current_data_sub, current_T_soil_C, U_csac, T_start_C);
%                 p = fminsearch(fun,max(current_data_sub.T_main_C));
% 
%                 T_prev_C = p(1);
% 
%                 for i=1:m
%                     T_prev_C = get_supply_temp(T_prev_C,Q_total(i),U_csac,x_diff_m(i),current_T_soil_C);
%                     if i==m
%                         T_main_fit_C(m) = T_prev_C;
%                     end
%                 end
%             end
%         end
%     end
% >>>>>>> parent of d5cf9ed (Improve temp models and filter stability)
