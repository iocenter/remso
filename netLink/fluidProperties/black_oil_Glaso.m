function [B_g,B_o,R_s] = black_oil_Glaso(p,R_sb,rho_g_sc,rho_o_sc,T)
% function [B_g,B_o,R_s] = black_oil_Glaso(p,R_sb,rho_g_sc,rho_o_sc,T)
%
% Computes the gas and oil formation volume factors B_g and B_o and the solution GOR R_s at a
% given pressure p and temperature T, bubble point GOR R_sb, and gas and oil densities rho_g_sc
% and rho_o_sc. The pressure p may be below or above the bubble point pressure.
%
% For the oil parameters p_b, B_o and R_s, use is made of the Glaso correlations, while for
% compressibility c_o and modified gas density rho_g_100, we used the Vazquez and Beggs
% correlations. 
%
% To compute the gas parameter B_g, use is made of the Sutton correlations for pseudo-critical
% pressure p_pc and temperature T_pc, and of the Dranchuk and Abu-Kassem approximation of the
% Standing-Katz correlation for the Z factor. 
% B_g = gas-formation volume factor, m^3/m^3
% B_o = oil-formation volume factor, m^3/m^3
% p = pressure, Pa 
% R_s = solution gas-oil ratio, m^3/m^3
% R_sb = solution gas-oil ratio at bubble point, m^3/m^3
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% T = temperature, deg. C
%
% JDJ, 08-01-02, last revised 10-05-13 

% Check input:
if R_sb > 363
    warning('R_sb above range of validity of Glaso black oil correlation.')
end
if T > 137
    warning('T above range of validity of Glaso black oil correlation.')
end

% Standard conditions:
p_sc = 100e3; % pressure at standard conditions, Pa
T_sc = 15; % temperature at standard conditions, deg. C

p_b = pres_bub_Glaso(R_sb,rho_g_sc,rho_o_sc,T); % bubble point pressure, Pa
if p < p_b % saturated oil 
    % Oil parameters:
    R_s = gas_oil_rat_Glaso(p,rho_g_sc,rho_o_sc,T);
    B_o = oil_form_vol_fact_Glaso(R_s,rho_g_sc,rho_o_sc,T);
    % Gas parameter:
    T_abs = T + 273.15; % absolute temperature, K
    p_pc = pres_pseu_crit_Sutton(rho_g_sc); % pseudo-critical pressure, Pa
    T_pc = temp_pseu_crit_Sutton(rho_g_sc); % pseudo-critical temperature, K
    p_pr = p / p_pc; % pseudo-reduced pressure, -
    T_pr = T_abs / T_pc;% pseudo reduced temperature, -
    Z = Z_factor_DAK(p_pr,T_pr); % Z factor, -
    B_g = gas_form_vol_fact(p,T_abs,Z);
else % undersaturated oil
    % Oil parameters:
    R_s = R_sb;
    B_ob = oil_form_vol_fact_Glaso(R_sb,rho_g_sc,rho_o_sc,T); % oil formation volume factor at
        % bubble point pressure, m^3/m^3
    rho_g_100 = rho_g_Vazquez_and_Beggs(p_sc,rho_g_sc,rho_o_sc,T_sc); % gas density at 100 psi,
        % kg/m^3
    c_o = compres_Vazquez_and_Beggs(p,R_sb,rho_g_100,rho_o_sc,T); % oil compressibiliy, 1/Pa
    B_o = oil_form_vol_fact_undersat(B_ob,c_o,p,p_b);
    % Gas parameter:
    B_g = 0; % No free gas.
end

