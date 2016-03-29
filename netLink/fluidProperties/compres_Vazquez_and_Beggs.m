function c_o = compres_Vazquez_and_Beggs(p,R_sb,rho_g_100,rho_o_sc,T)
% function c_o = compres_Vazquez_and_Beggs(p,R_sb,rho_g_100,rho_o_sc,T)
%
% Computes the compressibility with the Vazquez and Beggs correlation converted to SI units.
%
% c_o = oil compressibility, 1/Pa
% p = pressure, Pa 
% R_sb = solution gas oil ratio at bubble point pressure, m^3/m^3
% rho_g_100 = gas density at 100 psig, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% T = temperature, deg. C
%
% JDJ, 02-03-01, last revised 10-05-13
%{
Changes by Thiago and Codas
Make the function compatible with ADI objects


%}


help01 = 27.8 * R_sb;
help02 = 31 * T;
help03 = 959 * rho_g_100;
help04 = 1784000/rho_o_sc;
c_o = (-2541 + help01 + help02 - help03 + help04) ./ (1e5*p);
