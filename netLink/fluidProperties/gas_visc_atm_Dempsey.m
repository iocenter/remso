function mu_g_p_sc = gas_visc_atm_Dempsey(M,T)
% function mu_g_p_sc = gas_visc_atm_Dempsey(M,T)
%
% Calculates the gas viscosity at atmosperic pressure as a function of mollar mass M and
% temperature T in SI units.  
%
% Use is made of an expression of Dempsey (1965) to approximate the correlation of Carr,
% Kobayashi and Burrows (1954). 
%
% M = molar mass, kg/mol
% mu_g_p_sc = viscosity at atmospheric pressure, Pa s
% T = temperature, deg. C
%
% JDJ, 01-10-02, last revised 05-02-15

b0 =  1.16620808E-05;
b1 =  3.04342760E-08;
b2 =  6.84808007E-12;
b3 = -1.11626158E-04;
b4 = -1.25617746E-07;
b5 = -2.91397349E-10;
b6 =  4.64955375E-04;
b7 =  4.29044857E-07;
b8 =  1.28865249E-09;

mu_g_p_sc = b0 + b1*T + b2*T^2 + b3*M + b4*T*M + b5*T^2*M + b6*M^2 + b7*T*M^2 + b8*T^2*M^2;

