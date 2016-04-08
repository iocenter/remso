function mu_o = oil_visc_sat_B_and_R(mu_od,R_s)
% function mu_o = oil_visc_sat_B_and_R(mu_od,R_s)
%
% Computes the saturated-oil viscosity using the Beggs and Robinson correlation in SI units.
%
% mu_o = saturated-oil viscosity, Pa s
% mu_od = dead-oil viscosity, Pa s
% R_s = solution gas-oil ratio, m^3/m^3
%
% JDJ, 24-09-02, last revised 03-02-15
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%}
c = 5.44*(R_s./0.178+150).^-0.338;
mu_o = (10.715e-3*(R_s./0.178+100).^-0.515)*(mu_od*1e3).^c;