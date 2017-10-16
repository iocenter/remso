function mu_o = oil_viscosity(p,R_sb,rho_g_sc,rho_o_sc,T)
% function mu_o = oil_viscosity(p,R_sb,rho_g_sc,rho_o_sc,T)
%
% Computes the oil viscosity at given pressure, temperature, producing GOR and oil and gas
% densities at standard conditions. The pressure may be below or above the bubble point
% pressure. 
%
% For the dead-oil viscosity and the saturated-oil viscosity use is made of the Beggs and
% Robinson(1975) correlations, while for the undersaturated-oil viscosity we used the Vazquez
% and Beggs (1980) correlation. For the black oil properties we use the Standing (1952)
% correlations. 
%
% mu_o = oil viscosity, Pa s
% p = pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% R_sb = gas-oil ratio at bubble point pressure, m^3/m^3
% T = temperature, deg. C
%
% JDJ, 07-10-02, last revised 10-05-13
%
%{
%Changes by Codas: Vectorization
%Compatibility with ADI
%}

% Dead-oil viscosity:
mu_od = oil_visc_dead_B_and_R(rho_o_sc,T);

% Black oil properties:
p_b = pres_bub_Standing(R_sb,rho_g_sc,rho_o_sc,T); % bubble point pressure, Pa

mu_o = (p+R_sb)*0; % initialization

% Oil viscosity:
c1 = p<p_b;
if any(c1)
   R_s = gas_oil_rat_Standing(p(c1),rho_g_sc,rho_o_sc,T(c1)); % solution gas-oil ratio, m^3/m^3
   mu_o(c1) = oil_visc_sat_B_and_R(mu_od(c1),R_s); % saturated oil viscosity, Pa s
end
if (any(~c1))
   mu_ob = oil_visc_sat_B_and_R(mu_od(~c1),R_sb(~c1)); % oil viscosity at bubble point, Pa s
   mu_o(~c1) = oil_visc_undersat_V_and_B(mu_ob,p(~c1),p_b(~c1)); % undersaturated oil viscosity, Pa s
end
