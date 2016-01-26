function mu_g = gas_viscosity(p,rho_g_sc,T)
% viscosity mu_g = gas_viscosity(p,rho_g_sc,T)
%
% Calculates the gas viscosity as a function of pressure, temperature and gas density at
% standard conditions in SI units.  
%
% Use is made of the Dempsey (1965) approximations of the Carr, Kobayashi and Burrows (1954)
% correlations. 
%
% mu_g = gas viscosity, Pa s
% p = pressure, Pa
% rho_g_sc = gas density at standard condition, kg/m^3
% T = temperature, deg. C
%
% JDJ, 04-10-02, revised 08-04-15

M = from_kg_per_m3_to_molar_mass(rho_g_sc); % molar mass, kg/mol
mu_g_p_sc = gas_visc_atm_Dempsey(M,T); % gas viscosity at atmospheric pressure, Pa s
p_pc = pres_pseu_crit_Sutton(rho_g_sc); % pseudo-critical pressure, Pa
T_pc = temp_pseu_crit_Sutton(rho_g_sc); % pseudo-critical temperature, K
p_pr = p/p_pc; % pseudo-reduced pressure, - 
T_abs = T + 273.15; % absolute temperature, K
T_pr = T_abs/T_pc; % pseudo-reduced temperature, - 
f = gas_visc_ratio_Dempsey(p_pr,T_pr); % gas viscosity ratio, -
mu_g = f * mu_g_p_sc;

