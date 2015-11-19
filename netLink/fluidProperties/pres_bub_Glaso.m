function p_b = pres_bub_Glaso(R_sb,rho_g_sc,rho_o_sc,T)
% function p_b = pres_bub_Glaso(R_sb,rho_g_sc,rho_o_sc,T)
%
% Computes the bubble point pressure with the Glaso correlation.
%
% R_sb = solution gas-oil ratio at bubble point, m^3/m^3
% p_b = bubble point pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m3
% rho_o_sc = oil density at standard conditions, kg/m3
% T = temperature, deg. C
%
% JDJ, 24-03-03, last revised 10-05-13

% Convert parameters to field units:
gamma_g = from_kg_per_m3_to_gas_grav(rho_g_sc); % gas specific gravity wrt. air, -
gamma_API = from_kg_per_m3_to_deg_API(rho_o_sc); % API gravity, deg. API
R_sb_FU = from_m3_per_m3_to_ft3_per_bbl(R_sb); % solution gas-oil ratio at bubble pt., ft^3/bbl
T_FU = from_deg_C_to_deg_F(T); % temperature, deg. F

% Compute Glaso correlatation:
a = 0.816;
b = 0.172;
c = 0.989;
p_b_star = (R_sb_FU/gamma_g)^a * T_FU^b / gamma_API^c; % "correlating number"
p_b_FU = 10^(1.7669 + 1.7447*log10(p_b_star) - 0.30218*(log10(p_b_star))^2);
       % bubble point pressure, psi
       
% Convert to SI-units:
p_b = from_psi_to_Pa(p_b_FU);
