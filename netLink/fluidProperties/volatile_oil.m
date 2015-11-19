function [B_g,B_o,r_s,R_s] = volatile_oil(p,R_sb,rho_g_sc,rho_o_sc,T,vol_oil)
% function [B_g,B_o,r_s,R_s] = volatile_oil(p,R_sb,rho_g_sc,rho_o_sc,T,vol_oil)
%
% Computes the volatile oil parameters B_g, B_o, r_s and R_s at a given pressure p and
% temperature T, through two-dimensional interpolation (i.e. between p values and between T
% values) in a three-dimensional matrix. The pressure may be below or above the bubble point
% pressure.  
%
% vol_oil = three-dimensional array (matrix) with volatile oil parameters. The matrix contains
%           values of T, p, B_g, B_o, r_s, and R_s. The first dimension is used to loop over
%           the temperature values, i.e. it just contains integers 1, 2, …, m, where m is the
%           number of temperature values. The second dimension loops over the pressure values
%           (i.e. integers 1, 2, …, n, where n is the number of pressure values), and the
%           third dimension contains the values of T, p, and the four volatile oil parameters.         
% B_g = gas-formation volume factor, m^3/m^3
% B_o = oil-formation volume factor, m^3/m^3
% p = pressure, Pa 
% r_s = solution oil-gas ratio, m^3/m^3
% R_s = solution gas-oil ratio, m^3/m^3
% R_sb = solution gas-oil ratio at bubble point pressure, m^3/m^3
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3

% T = temperature, deg. C
%
% JDJ, 08-05-13, last revised 27-10-13

% Check input:
lb = 0.999; % lower bound (multiplier), -
ub = 1.001; % upper bound (multiplier), -
[m,n,k] = size(vol_oil);
T_lo = vol_oil(1,1,1);
p_lo = vol_oil(1,1,2);
T_hi = vol_oil(m,1,1);
p_hi = vol_oil(1,n,2);
if T < T_lo
    T
    error('T below range of validity of volatile oil table.')
end
if T > T_hi
    T
    error('T above range of validity of volatile oil table.')
end
if p < p_lo
    p
    error('p below range of validity of volatile oil table.')
end
if p > p_hi
    p
    error('p above range of validity of volatile oil table.')
end
if rho_g_sc < lb*0.80 || rho_g_sc > ub*0.80
    rho_g_sc
    error('Value of rho_g_sc is inconsistent with value used to generate volatile oil table (0.80 kg/m^3).')
end
if rho_o_sc < lb*800 || rho_o_sc > ub*800 
    rho_o_sc
    error('Value of rho_o_sc is inconsistent with value used to generate volatile oil table (800 kg/m^3).')
end
if R_sb < lb*450 || R_sb > ub*450  
    R_sb
    error('Value of R_sb is inconsistent with value used to generate volatile oil table (450 m^3/m^3).')
end

% Determine the volatile oil parameters through interpolation:
B_g = interp2(vol_oil(:,:,1)',vol_oil(:,:,2)',vol_oil(:,:,3)',T,p);
B_o = interp2(vol_oil(:,:,1)',vol_oil(:,:,2)',vol_oil(:,:,4)',T,p);
r_s = interp2(vol_oil(:,:,1)',vol_oil(:,:,2)',vol_oil(:,:,5)',T,p);
R_s = interp2(vol_oil(:,:,1)',vol_oil(:,:,2)',vol_oil(:,:,6)',T,p);
