function p_b = pres_bub_Standing(R_sb,rho_g_sc,rho_o_sc,T)
% function p_b = pres_bub_Standing(R_sb,rho_g_sc,rho_o_sc,T)
%
% Computes the bubble point pressure with a Standing correlation converted to SI units.
%
% R_sb = gas-oil ratio at bubble point pressure, m^3/m^3
% p_b = bubble point pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m3
% rho_o_sc = oil density at standard conditions, kg/m3
% T = temperature, deg. C
%
% JDJ, 27-02-01, last revised 10-05-13
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%}
% check for presence of gas:
if rho_g_sc == 0
    p_b = 1.e5; % atmospheric pressure
else
    help01 = (10.^(0.00164*T))/(10.^(1768./rho_o_sc));
    p_b = 125e3 * ((716*R_sb./rho_g_sc).^0.83 * help01 - 1.4);
end

% Reality check:
if p_b < 1.e5
    p_b = 1.e5; % atmospheric pressure
end
