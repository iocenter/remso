function B_o = oil_form_vol_fact_Standing(R_s,rho_g_sc,rho_o_sc,T)
% function B_o = oil_form_vol_fact_Standing(R_s,rho_g_sc,rho_o_sc,T)
%
% Computes the oil formation volume factor with a Standing correlation converted to SI units.
%
% B_o = oil formation volume factor, m^3/m^3
% R_s = solution gas-oil ratio, m^3/m^3
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% T = temperature, deg. C
%
% JDJ, 27-02-01, last revised 10-05-13
%{
Changes by Thiago and Codas
Make the function compatible with ADI objects
%}
help01 = sqrt(rho_g_sc./rho_o_sc); 
B_o = 0.9759 + 12e-5 *(160 * R_s * help01 + 2.25 * T + 40).^1.2;