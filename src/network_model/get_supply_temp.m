%%% getSupplyTemp
% finds the expected supply temperature based on temperature in the main
% pipe and knowledge on flow, insulation, ambient temperature and length of
% the pipe.
% The function below assumes that there is a linear temperature profile in
% the subpipe. As such, it only works, when ther is sufficient flow.

function T_supply = get_supply_temp(T_main_C, flow_kg_h, U,length_pipe_m,T_soil_C)

    c = 4.186e3; % J/kg/°C

    flow_kg_s = reshape(flow_kg_h,[],1)./3600; % translate from kg^3/hour to kg^3/s (SI-units)
    U = U(:); % W/K/m
    length_pipe_m = length_pipe_m(:); % m

    T_main_C = T_main_C(:);

    B = U.*length_pipe_m; % W/K
    theta = B/(2*c);
    
    alpha = (flow_kg_s - theta) ./ (flow_kg_s + theta);

    term1 = alpha .* T_main_C;
    term2 = (1-alpha) .* T_soil_C;

    T_supply = term1+term2;