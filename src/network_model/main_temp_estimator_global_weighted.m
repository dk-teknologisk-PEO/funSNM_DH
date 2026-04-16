function cost = main_temp_estimator_global_weighted(params, current_data, valid_idx, weights, T_soil_C, U_csac)
    T_start = params(1);
    master_offset = params(2);
    
    Q_total = cumsum(current_data.flow_kg_h, 'reverse');
    x_diff_m = current_data.x_pos_m;
    x_diff_m(2:end) = diff(x_diff_m);
    
    T_prev = T_start;
    T_main_sim = zeros(size(current_data,1), 1);
    for i = 1:size(current_data,1)
        T_prev = get_supply_temp(T_prev, Q_total(i), U_csac, x_diff_m(i), T_soil_C);
        T_main_sim(i) = T_prev;
    end
    
    % Weighted residuals
    residuals = T_main_sim(valid_idx) - current_data.T_main_C(valid_idx) - master_offset;
    cost = sum(weights .* residuals.^2);
end