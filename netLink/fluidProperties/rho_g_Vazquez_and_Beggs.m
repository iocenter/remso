function rho_g_100 = rho_g_Vazquez_and_Beggs(p_sep,rho_g_sep,rho_o_sc,T_sep)
% rho_g_100 = rho_g_Vazquez_and_Beggs(p_sep,rho_g_sep,rho_o_sc,T_sep)
%
% Computes the equivalent gas density as if determined from a sample taken at a separator
% pressure of 689 kPa (100 psi). Input is the gas density rho_g determined from a sample
% taken at another (separator) pressure p_sep and temperature T_sep. Use is made of a
% correlation from Vazquez and Beggs, converted to SI units. 
%
% p_sep = separator pressure, Pa
% rho_g_sep = gas density at p_sep, kg/m^3  
% rho_g_100 = gas density at 100 psi, kg/m^3  
% rho_o_sc = oil density at standard conditions, kg/m^3
% T_sep = separator temperature, deg. C
%
% JDJ, 02-03-01, last revised 05-02-15

help01 = 141500 / rho_o_sc - 131.5;
help02 = 1.8 * T_sep + 32;
help03 = p_sep /790.8e3; 
rho_g_100 = rho_g_sep * (1 + 5.912e-5 * help01 * help02 * log10(help03));