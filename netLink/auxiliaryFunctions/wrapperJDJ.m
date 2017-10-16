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
rho_sc = [streams(1).gas_dens, streams(1).oil_dens, streams(1).water_dens]; % in kg/m^3 

assert(all(vertcat(streams.gas_dens)==streams(1).gas_dens));
assert(all(vertcat(streams.oil_dens)==streams(1).oil_dens));
assert(all(vertcat(streams.water_dens)==streams(1).water_dens));

%%TODO: consider fluid properties (stream) as a feature of the network or vectorize rho_sc

s_in = 0; 
s_out = 1; 
T = vertcat(pipelines.temp) - 273.15; 
T_in = T;
T_out =T;
