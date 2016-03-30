function R_s = gas_oil_rat_Standing(p,rho_g_sc,rho_o_sc,T)
% function R_s = gas_oil_rat_Standing(p,rho_g_sc,rho_o_sc,T)
%
% Computes the solution gas-oil ratio with a Standing correlation converted to SI units. 
% 
% R_s = solution gas-oil ratio, m^3/m^3
% p = pressure, Pa
% rho_g_sc, gas density at standard conditions, kg/m^3  
% rho_o_sc, oil density at standard conditions, kg/m^3  
% T = temperature, deg. C
%
% JDJ, 27-02-01, last revised 10-05-13
%{
Changes Codas
Make the function compatible with ADI objects
%}
help01 = 10.^(1768./rho_o_sc - 0.00164*T);
R_s = (rho_g_sc./716)*((8e-6*p+1.4)*help01).^1.2048;