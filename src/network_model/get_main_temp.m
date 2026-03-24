    %%% getMainTemp.m %%%
% function, taking in an array of meter temperatures and -flow rates and
% finding the temperature in the main pipe.
% The temperatures must be in °C and the flow in ltr/hour. The meter
% measurement is expected to be combined in a matrix with temperatures in
% the first column and corresponding flows in the second column.
% U-value of the subpipe is given as an array, as well as the length of the
% subpipes, L.
% A common ambient temperature is expected for all measurements.


function T_main = get_main_temp(T_supply_C, flow_kg_h,U ,L_pipe_m ,T_soil_C)
    
    c = 4.186e3; % j/kg/K

    flow_kg_s = flow_kg_h./3600;
    B = U.*L_pipe_m;
    nom = T_supply_C.*(c*flow_kg_s + (B./2)) - B.*T_soil_C;
    denom = c*flow_kg_s-(B./2);
    T_main = nom./denom;