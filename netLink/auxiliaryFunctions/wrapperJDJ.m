function [ s,alpha,d,e,oil, rho_sc,s_in,s_out,T, T_in,T_out] = wrapperJDJ(E)
%wrapperJDJ Wrapper for JDJ dp function
s = 0;
pipelines = vertcat(E.pipeline);
streams = vertcat(E.stream);
if any((vertcat(pipelines.ang) > pi/2)) || any((vertcat(pipelines.ang) < -pi/2))
    warning('check angle definition')
end
alpha = pi/2-vertcat(pipelines.ang);
d = vertcat(pipelines.diam);
e = vertcat(pipelines.roughness);
oil = 1; % black oil; parameters computed with the aid of Standing correlations
rho_sc = [streams.gas_dens, streams.oil_dens, streams.water_dens]; % in kg/m^3 
s_in = 0; 
s_out = 1; 
T = vertcat(pipelines.temp) - 273.15; 
T_in = T;
T_out =T;
