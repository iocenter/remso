function mu_o = oil_visc_undersat_V_and_B(mu_ob,p,p_b)
% function mu_o = oil_visc_undersat_V_and_B(mu_ob,p,p_b)
%
% Computes the undersaturated-oil viscosity using the Vazquez and Beggs correlation in SI
% units. 
%
% mu_o = undersaturated-oil viscosity, Pa s
% mu_ob = oil viscosity at bubble point, Pa s
% p = pressure, Pa
% p_b = bubble point pressure, Pa
%
% JDJ, 27-09-02, last revised 10-05-13

d = 7.2e-5*p^1.187*exp(-11.513-1.30e-8*p);
mu_o = mu_ob*(p/p_b)^d;