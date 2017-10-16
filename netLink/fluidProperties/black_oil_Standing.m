function [B_g,B_o,R_s,Z] = black_oil_Standing(p,R_sb,rho_g_sc,rho_o_sc,T,hasSurfaceGas,Z)
% function [B_g,B_o,R_s] = black_oil_Standing(p,R_sb,rho_g_sc,rho_o_sc,T)
%
% Computes the gas and oil formation volume factors B_g and B_o and the solution GOR R_s at a
% given pressure p and temperature T, bubble point GOR R_sb, and gas and oil densities rho_g_sc
% and rho_o_sc. The pressure p may be below or above the bubble point pressure.
%
% For the oil parameters p_b, B_o and R_s, use is made of the Standing correlations, while for
% compressibility c_o and modified gas density rho_g_100, we used the Vazquez and Beggs
% correlations. 
%
% To compute the gas parameter B_g, use is made of the Sutton correlations for pseudo-critical
% pressure p_pc and temperature T_pc, and of the Dranchuk and Abu-Kassem approximation of the
% Standing-Katz correlation for the Z factor.
%
% B_g = gas-formation volume factor, m^3/m^3
% B_o = oil-formation volume factor, m^3/m^3
% p = pressure, Pa 
% R_sb = solution gas-oil ratio at bubble point pressure, m^3/m^3
% R_s = solution gas-oil ratio, m^3/m^3
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% T = temperature, deg. C
%
% JDJ, 08-01-02, last revised 09-05-13 

if nargin < 6
    hasSurfaceGas = true;
end
if nargin < 7
    Z = [];
end

% Check input:
if any(R_sb > 254)
    warning('R_sb above range of validity of Standing black oil correlation.')
end
if any(T > 125)
    warning('T above range of validity of Standing black oil correlation.')
end
if any(rho_o_sc < 725)
    warning('rho_o_sc below range of validity of Standing black oil correlation.')
end
if any(rho_o_sc > 956)
    warning('rho_o_sc above range of validity of Standing black oil correlation.')
end
if any(rho_g_sc < 0.73)
    warning('rho_g_sc below range of validity of Standing black oil correlation.')
end
if any(rho_g_sc > 1.17)
    warning('rho_g_sc above range of validity of Standing black oil correlation.')
end

% Standard conditions:
p_sc = 100e3; % pressure at standard conditions, Pa
T_sc = 15; % temperature at standard conditions, deg. C

p_b = pres_bub_Standing(R_sb,rho_g_sc,rho_o_sc,T); % bubble point pressure, Pa

R_s = (p+R_sb)*0; %% initialization
B_o = R_s;  %% initialization 

% Oil parameters:
c1 = p <= p_b;
if any(c1) % saturated oil 
    R_s(c1) = gas_oil_rat_Standing(p(c1),rho_g_sc,rho_o_sc,T(c1));
    B_o(c1) = oil_form_vol_fact_Standing(R_s(c1),rho_g_sc,rho_o_sc,T(c1));
end
if any(~c1) % undersaturated oil
    R_s(~c1) = R_sb(~c1);
    B_ob = oil_form_vol_fact_Standing(R_s(~c1),rho_g_sc,rho_o_sc,T(~c1)); % oil formation volume factor
    % at bubble point pressure, m^3/m^3
    rho_g_100 = rho_g_Vazquez_and_Beggs(p_sc,rho_g_sc,rho_o_sc,T_sc); % gas density at 100 psi,
    % kg/m^3
    c_o = compres_Vazquez_and_Beggs(p(~c1),R_s(~c1),rho_g_100,rho_o_sc,T(~c1)); % oil compressibility, 1/Pa
    B_o(~c1) = oil_form_vol_fact_undersat(B_ob,c_o,p(~c1),p_b(~c1));
end
% Gas parameter:
B_g = B_o; %initiallization
c1 = hasSurfaceGas | R_s > 0;
if any(c1)
    T_abs = T(c1) + 273.15; % absolute temperature, K
    p_pc = pres_pseu_crit_Sutton(rho_g_sc); % pseudo-critical pressure, Pa
    T_pc = temp_pseu_crit_Sutton(rho_g_sc); % pseudo-critical temperature, K
    p_pr = p(c1) / p_pc; % pseudo-reduced pressure, -
    T_pr = T_abs / T_pc; % pseudo reduced temperature, -
    if nargin < 7 || isempty(Z) || any(isnan(Z))
        Z = Z_factor_DAK(p_pr,T_pr); % Z factor, -
    else
        Z = Z_factor_DAK(p_pr,T_pr,Z); % for the gradient correction!
    end
    B_g(c1) = gas_form_vol_fact(p(c1),T_abs,Z);
end
if any(~c1)
    B_g(~c1) = nan;    
end
