function [dpds_tot,f_work,Z] =  dpBeggsBrillJDJ(E, qoE, qwE, qgE, pV,hasSurfaceGas,f_work,Z)
%% dpBeggsBrillJDJ - wrapper to Beggs_Brill_dpds
%% dp: calculates the pressure drop for the average pressure pV
%%     the functions is implemented for SI units (input).

%TODO: check angles and flow direction conventions.
%TODO: interface stream parameters such as viscosity and specific gas
%gravity. Check E.stream

s = 0;
if any((E.pipeline.ang > pi/2)) || any((E.pipeline.ang < -pi/2))
    warning('check angle definition')
end
if nargin < 7
    f_work = [];
    Z = [];
end
if nargin < 6
    hasSurfaceGas = true;
end
alpha = pi/2-E.pipeline.ang;
d = E.pipeline.diam;
e = E.pipeline.roughness;
oil = 1; % black oil; parameters computed with the aid of Standing correlations
rho_sc = [E.stream.gas_dens,E.stream.oil_dens,E.stream.water_dens]; % in kg/m^3 
s_in = 0; 
s_out = 1; 
T = E.pipeline.temp - 273.15; 
T_in = T;
T_out =T;

[dpVector,f_work,Z] = Beggs_Brill_dpds(s,pV,[],alpha,d,e,oil,qgE, qoE, qwE,       rho_sc,s_in,s_out,T_in,T_out,hasSurfaceGas,f_work,Z);
dpds_tot = dpVector(1);
end
