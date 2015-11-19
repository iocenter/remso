function R_s = gas_oil_rat_Glaso(p,rho_g_sc,rho_o_sc,T)
% function R_s = gas_oil_rat_Glaso(p,rho_g_sc,rho_o_sc,T)
%
% Computes the solution gas-oil ratio with the Glaso correlation. This correlation has been
% obtained through inversion of the Glaso bubble point pressure correlation.
%
% R_s = solution gas-oil ratio, m^3/m^3
% p = pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m3
% rho_o_sc = oil density at standard conditions, kg/m3
% T = temperature, deg. C
%
% JDJ, 24-03-03, last revised 10-05-13

% Convert parameters to field units:
gamma_g = from_kg_per_m3_to_gas_grav(rho_g_sc); % gas specific gravity wrt. air, -
gamma_API = from_kg_per_m3_to_deg_API(rho_o_sc); % API gravity, deg. API
p_FU = from_Pa_to_psi(p); % pressure, psi
T_FU = from_deg_C_to_deg_F(T); % temperature, deg. F

% Compute correlation:
p_b_star = 10^(2.8869-(14.1811-3.3039*log10(p_FU))^0.5);
R_s_FU = gamma_g*(p_b_star*gamma_API^0.989/T_FU^0.172)^1.2255; % ft^3/bbl

% Convert to SI-units:
R_s = from_ft3_per_bbl_to_m3_per_m3(R_s_FU);