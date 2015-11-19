function B_o = oil_form_vol_fact_Glaso(R_s,rho_g_sc,rho_o_sc,T)
% function B_o = oil_form_vol_fact_Glaso(R_s,rho_g_sc,rho_o_sc,T)
%
% Computes the oil formation volume factor with the Glaso correlation.
%
% B_o = oil formation volume factor, m^3/m^3
% R_s = solution gas-oil ratio, m^3/m^3
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% T = temperature, deg. C
%
% JDJ, 24-03-03, last revised 10-05-13

% Convert parameters to field units:
gamma_g = from_kg_per_m3_to_gas_grav(rho_g_sc); % gas specific gravity w.r.t. air, -
gamma_o = from_kg_per_m3_to_liq_grav(rho_o_sc); % oil specific gravity w.r.t. water, -
R_s_FU = from_m3_per_m3_to_ft3_per_bbl(R_s);
T_FU = from_deg_C_to_deg_F(T);

% Compute Glaso correlatation:
B_ob_star = R_s_FU*(gamma_g/gamma_o)^0.526 + 0.968*T_FU; % "correlating number"
A = -6.58511 + 2.91329*log10(B_ob_star) - 0.27683*(log10(B_ob_star))^2;
B_o_FU = 1 + 10^A;

% Convert to SI-units:
B_o = B_o_FU; % trivial
