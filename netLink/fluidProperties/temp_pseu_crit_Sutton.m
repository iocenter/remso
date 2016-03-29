function T_pc = temp_pseu_crit_Sutton(rho_g_sc)
% function T_pc = temp_pseu_crit_Sutton(rho_g_sc)
%
% Calculates the pseudo-critical temperature of a gas mixture
% with unknown composition, using the Sutton (1985) correlation
% converted to SI units.
% 
% rho_g_sc, gas density at standard conditions, kg/m^3
% T_pc = pseudo-critical temperature, K
%
% JDJ, 05-03-01, last revised 14-04-12
%
%{
Changes by Thiago and Codas
Make the function compatible with ADI objects
%}
T_pc = 94.0 + 157.9 * rho_g_sc - 27.2 * rho_g_sc.^2;