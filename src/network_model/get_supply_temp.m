%%% getSupplyTemp
% finds the expected supply temperature based on temperature in the main
% pipe and knowledge on flow, insulation, ambient temperature and length of
% the pipe.
% The function below assumes that there is a linear temperature profile in
% the subpipe. As such, it only works, when ther is sufficient flow.

function T_supply = get_supply_temp(T_main_C, flow_kg_h, U,length_pipe_m,T_soil_C)
    
    flow_kg_s = reshape(flow_kg_h,[],1)./3600; % translate from kg^3/hour to kg^3/s (SI-units)
    U = reshape(U,[],1); % W/K/m
    length_pipe_m = reshape(length_pipe_m,[],1); % m
    c = 4.186e3; % J/kg/°C
    B = U.*length_pipe_m; % W/K
    theta = B/(2*c);
    
    term1 = ((flow_kg_s-theta)./(flow_kg_s+theta)).*T_main_C;
    term2 = ((2*theta)./(flow_kg_s+theta)).*T_soil_C;

    T_supply = term1+term2;