function csac_state = apply_master_offset(csac_state, ukf_master_offset, pf_master_offset, ...
        num_active_houses, params)
%APPLY_MASTER_OFFSET Apply damped master offset to all houses.

    if ~isnan(ukf_master_offset) && abs(ukf_master_offset) > params.master_offset_deadzone
        if num_active_houses >= params.min_active_houses
            damped = params.master_offset_gamma * ukf_master_offset;
            for i = 1:csac_state.num_houses
                csac_state.ukf_states{i}.x(1) = csac_state.ukf_states{i}.x(1) - damped;
            end
        end
    end

    if ~isnan(pf_master_offset) && abs(pf_master_offset) > params.master_offset_deadzone
        if num_active_houses >= params.min_active_houses
            damped = params.master_offset_gamma * pf_master_offset;
            for i = 1:csac_state.num_houses
                csac_state.pf_states{i}.x(1) = csac_state.pf_states{i}.x(1) - damped;
                csac_state.pf_particles{i}(1,:) = csac_state.pf_particles{i}(1,:) - damped;
            end
        end
    end

    if params.error_meaning && ~params.debug_disable_mean_centering
        mean_ukf = mean(csac_state.ukf_offsets, 'omitmissing');
        mean_pf = mean(csac_state.pf_offsets, 'omitmissing');
        for i = 1:csac_state.num_houses
            if ~isnan(mean_ukf) && ~isnan(csac_state.ukf_states{i}.x(1))
                csac_state.ukf_states{i}.x(1) = csac_state.ukf_states{i}.x(1) - mean_ukf;
            end
            if ~isnan(mean_pf) && ~isnan(csac_state.pf_states{i}.x(1))
                csac_state.pf_states{i}.x(1) = csac_state.pf_states{i}.x(1) - mean_pf;
                csac_state.pf_particles{i}(1,:) = csac_state.pf_particles{i}(1,:) - mean_pf;
            end
        end
    end
end