function p_pc = pres_pseu_crit_Sutton(rho_g_sc)
% function p_pc = pres_pseu_crit_Sutton(rho_g_sc)
%
% Calculates the pseudo-critical pressure of a gas mixture
% with unknown composition, using the Sutton (1985) correlation
% converted to SI units.
% 
% p_pc = pseudo-critical pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m^3
%
% JDJ, 05-03-01, last revised 13-04-12
%{
Changes by Thiago and Codas
Make the function compatible with ADI objects
%}

p_pc = 5218e3 - 734e3 * rho_g_sc - 16.4e3 * rho_g_sc.^2;