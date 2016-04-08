function B_gw = gas_form_vol_fact(p,T_abs,Z)
% function B_gw = gas_form_vol_fact(p,T_abs,Z)
%
% Computes the wet gas formation volume factor in SI units.
%
% B_gw = wet gas formation volume factor, m^3/m^3
% p = presssure, Pa
% T = temperature, K
% Z = gas compressibility factor, -
%
% JDJ, 02-01-02, last revised 10-05-13 
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%}

p_sc = 100e3; % pressure at standard conditions, Pa
T_sc_abs = 15 + 273.15; % temperature at standard conditions, K
Z_sc = 1; % gas compressibility factor at standard conditions, -
B_gw = (p_sc * T_abs * Z) ./ (p * T_sc_abs * Z_sc);
