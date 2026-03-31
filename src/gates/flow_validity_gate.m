function [is_valid, alpha, theta, min_flow_kg_h] = flow_validity_gate(flow_kg_h, U, length_pipe_m, alpha_min)
%FLOW_VALIDITY_GATE Checks if flow is sufficiently high for model validity.
%
% Inputs:
%   flow_kg_h      - flow in kg/h
%   U              - pipe U-value [W/(m K)]
%   length_pipe_m  - pipe length [m]
%   alpha_min      - minimum allowed (m_dot-theta)/(m_dot+theta)
%
% Outputs:
%   is_valid       - true if the criterion is satisfied
%   alpha          - current value of (m_dot-theta)/(m_dot+theta)
%   theta          - thermal parameter theta
%   min_flow_kg_h  - minimum flow corresponding to alpha_min

    c = 4186; % J/(kg K)
    flow_kg_s = flow_kg_h / 3600;

    theta = (U * length_pipe_m) / (2 * c);

    % avoid divide-by-zero pathology
    if (flow_kg_s + theta) <= 0
        alpha = -Inf;
        is_valid = false;
        min_flow_kg_h = Inf;
        return
    end

    alpha = (flow_kg_s - theta) / (flow_kg_s + theta);

    % exact minimum flow from alpha > alpha_min
    % m_dot > theta * (1 + alpha_min)/(1 - alpha_min)
    min_flow_kg_s = theta * (1 + alpha_min) / (1 - alpha_min);
    min_flow_kg_h = min_flow_kg_s * 3600;

    is_valid = alpha > alpha_min;
end