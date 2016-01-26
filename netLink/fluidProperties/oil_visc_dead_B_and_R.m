function mu_od = oil_visc_dead_B_and_R(rho_o_sc,T)
% function mu_od = oil_visc_dead_B_and_R(rho_o_sc,T)
%
% Computes the dead-oil viscosity using the Beggs and Robinson correlation
% in SI units.
%
% mu_od = dead-oil viscosity, Pa s
% rho_o_sc = oil density at standard conditions, kg/m^3
% T = temperature, C
%
% JDJ, 24-09-02, last revised 10-05-13

b = 5.693-2.863*10^3/rho_o_sc;
a = 10^b / (1.8*T+32)^1.163;
mu_od = 10^-3*(10^a-1);