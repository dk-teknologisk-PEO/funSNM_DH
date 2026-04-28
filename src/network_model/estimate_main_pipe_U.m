function [U_main_new, diagnostics] = estimate_main_pipe_U(all_cs, csac_ids, topology, U_main_current, config)
%ESTIMATE_MAIN_PIPE_U Estimates main pipe U-value from position-dependent
%   patterns in per-CSAC average offsets.
%
%   If U_main is wrong, the average offset per CSAC will show a systematic
%   trend with CSAC position on the main pipe.

    cfg = config.project.main_pipe_U_estimation;
    num_csacs = numel(csac_ids);

    %% Extract per-CSAC average offsets and positions
    avg_offsets = nan(num_csacs, 1);
    positions = zeros(num_csacs, 1);
    num_converged = zeros(num_csacs, 1);

    P_thresh = cfg.P_offset_convergence_threshold^2;

    for c = 1:num_csacs
        cs = all_cs{c};

        % CSAC position on main pipe
        csac_topo_idx = find([topology.cul_de_sacs.id] == csac_ids(c));
        if ~isempty(csac_topo_idx)
            positions(c) = topology.cul_de_sacs(csac_topo_idx).dist_on_main_m;
        end

        % Compute mean offset for converged houses
        offsets = zeros(cs.num_houses, 1);
        valid = false(cs.num_houses, 1);
        for i = 1:cs.num_houses
            offsets(i) = cs.ukf_states{i}.x(1);
            if cs.ukf_states{i}.P(1,1) < P_thresh && cs.ukf_states{i}.P(1,1) > 0
                valid(i) = true;
            end
        end

        if sum(valid) >= 2
            avg_offsets(c) = mean(offsets(valid));
            num_converged(c) = sum(valid);
        end
    end

    %% Check if enough CSACs have valid data
    valid_csacs = isfinite(avg_offsets);
    if sum(valid_csacs) < 2
        U_main_new = U_main_current;
        diagnostics.adjusted = false;
        diagnostics.reason = 'insufficient_csacs';
        return;
    end

    %% Weighted linear regression: avg_offset vs position
    w = zeros(num_csacs, 1);
    w(valid_csacs) = num_converged(valid_csacs);
    w = w / sum(w);

    % De-mean offsets (remove network-wide shared offset)
    avg_offsets_dm = avg_offsets - sum(w .* avg_offsets);

    [slope, correlation] = weighted_linear_slope(positions, avg_offsets_dm, w);

    pos_range = max(positions) - min(positions);
    total_gradient = slope * pos_range;

    %% Adjust U_main
    if abs(correlation) < cfg.min_correlation
        U_main_new = U_main_current;
        diagnostics.adjusted = false;
        diagnostics.reason = 'weak_correlation';
        diagnostics.slope = slope;
        diagnostics.correlation = correlation;
        diagnostics.total_gradient = total_gradient;
        return;
    end

    % Positive gradient → U_main too low → increase
    raw_adjustment = cfg.gain * total_gradient;
    adjustment = max(-cfg.max_adjustment_per_step, min(cfg.max_adjustment_per_step, raw_adjustment));

    U_main_new = U_main_current + adjustment;
    U_main_new = max(cfg.U_min, min(cfg.U_max, U_main_new));

    %% Diagnostics
    diagnostics.adjusted = true;
    diagnostics.reason = 'updated';
    diagnostics.slope = slope;
    diagnostics.correlation = correlation;
    diagnostics.total_gradient = total_gradient;
    diagnostics.adjustment = adjustment;
    diagnostics.avg_offsets = avg_offsets;
    diagnostics.positions = positions;
    diagnostics.U_old = U_main_current;
    diagnostics.U_new = U_main_new;
end


function [slope, correlation] = weighted_linear_slope(x, y, w)
    x_mean = sum(w .* x);
    y_mean = sum(w .* y);
    x_c = x - x_mean;
    y_c = y - y_mean;
    num = sum(w .* x_c .* y_c);
    den = sum(w .* x_c.^2);
    if abs(den) < 1e-12
        slope = 0;
        correlation = 0;
        return;
    end
    slope = num / den;
    std_x = sqrt(sum(w .* x_c.^2));
    std_y = sqrt(sum(w .* y_c.^2));
    if std_x > 0 && std_y > 0
        correlation = num / (std_x * std_y);
    else
        correlation = 0;
    end
end