function B_gw = gas_form_vol_fact_direct(p,rho_g_sc,T_abs)
% function B_gw = gas_form_vol_fact_direct(p,rho_g_sc,T_abs)
%
% Computes the wet gas formation volume factor in SI units. Equal to gas_form_vol_fact.m but
% takes different arguments to avoid the need to compute intermediate steps (pseudo properties
% and Z factor).  
%
% B_gw = wet gas formation volume factor, m^3/m^3
% p = pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m^3
% T_abs = absolute temperature, K
%
% JDJ, 27-10-13, last revised 27-10-13 

Z = Z_factor_DAK_direct(p,rho_g_sc,T_abs); % Z factor, -
p_sc = 100e3; % pressure at standard conditions, Pa
T_sc_abs = 15 + 273.15; % temperature at standard conditions, K
Z_sc = 1; % gas compressibility factor at standard conditions, -
B_gw = (p_sc * T_abs * Z) / (p * T_sc_abs * Z_sc);
