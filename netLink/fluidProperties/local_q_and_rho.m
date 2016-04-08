function [q_g,q_o,q_w,rho_g,rho_o,rho_w,Z] = local_q_and_rho(oil,p,q_g_sc,q_o_sc,q_w_sc,R_sb,rho_sc,T,hasSurfaceGas,Z)
% function [q,rho] = local_q_and_rho(oil,p,q_sc,R_sb,rho_sc,T)
%
% Computes the local values of q = [q_g,q_o,q_w] and rho = [rho_g,rho_o,rho_w] from
% q_sc = [q_g_sc,q_o_sc,q_w_sc] and rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc] at a given
% pressure p and temperature T. A choice can be made between using a black oil model or a
% volatile oil table with the aid of parameter 'oil'. 
%
% oil = parameter to select black oil model or volatile oil table, -  
%   oil = 1: black oil; parameters computed with the aid of Standing correlations
%   oil = 2: black oil; parameters computed with the aid of Glaso correlations
%   oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
% p = pressure, Pa
% q = [q_g,q_o,q_w]
% q_g = gas flow rate at local conditions, m^3/s
% q_o = oil flow rate at local conditions, m^3/s
% q_w = water flow rate at local conditions, m^3/s
% q_sc = [q_g_sc,q_o_sc,q_w_sc]
% q_g_sc = gas flow rate at standard conditions, m^3/s
% q_o_sc = oil flow rate at standard conditions, m^3/s
% q_w_sc = water flow rate at standard conditions, m^3/s
% R_sb = gas-oil ratio at bubble point pressure, m^3/m^3
% rho = [rho_g,rho_o,rho_w]
% rho_g = gas density at local conditions, kg/m^3
% rho_o = oil density at local conditions, kg/m^3
% rho_w = water density at local conditions, kg/m^3
% rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc]
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% rho_w_sc = water density at standard conditions, kg/m^3
% T = temperature, deg. C
%
% JDJ, 05-02-02, last revised 09-05-13
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%vectorization
%}
% Compute auxiliary variables:

if nargin < 9
    hasSurfaceGas = true;
end
if nargin < 10
    Z = [];
end

rho_g_sc = rho_sc(1);
rho_o_sc = rho_sc(2);

% Compute oil parameters:
%   B_g = gas-formation volume factor, m^3/m^3
%   B_o = oil-formation volume factor, m^3/m^3
%   r_s = solution oil-gas ratio, m^3/m^3 (not relevant for black oil)
%   R_s = solution gas-oil ratio, m^3/m^3
switch oil
    case 1
        [B_g,B_o,R_s,Z] = black_oil_Standing(p,R_sb,rho_g_sc,rho_o_sc,T,hasSurfaceGas,Z);
        r_s = 0;
    case 2
        [B_g,B_o,R_s] = black_oil_Glaso(p,R_sb,rho_g_sc,rho_o_sc,T);
        r_s = 0;
    case 3
        %file_name = 'vol_oil_table_01';
        %load(file_name,'vol_oil');
        %[B_g,B_o,r_s,R_s] = volatile_oil(p,R_sb,rho_g_sc,rho_o_sc,T,vol_oil);
        B_g = 0;
        B_o = 0;
        r_s = 0;
        R_s = 0;
        assert(hasSurfaceGas);
end

% Assemble transformation matrices T_q and T_rho: 
D = 1 - R_s*r_s;

 computeLocalGas = hasSurfaceGas | R_s > 0;
 if any(computeLocalGas)
    zeroVal = (p+q_g_sc+q_o_sc+q_w_sc)*0;
    q_g = zeroVal;
    rho_g = zeroVal;
 else
    q_g = zeros(size(double(q_g_sc)));
    rho_g = zeros(size(double(q_g_sc)));
 end

 %T_q = cell(3,3);
% if computeLocalGas
% T_q{1,1} =  B_g./D;
% T_q{1,2} = -B_g*R_s./D;
% T_q{1,3} =  0;
% end
% T_q{2,1} = -B_o*r_s./D;
% T_q{2,2} =  B_o./D;
% T_q{2,3} =  0;
% T_q{3,1} =  0;
% T_q{3,2} =  0;
% T_q{3,3} =  1;
% 
% T_rho = cell(3,3);
% 
% if computeLocalGas
%     T_rho{1,1} = 1./B_g;
%     T_rho{1,2} = r_s./B_g;
%     T_rho{1,3} = 0;
% end
% T_rho{2,1} = R_s./B_o;
% T_rho{2,2} = 1./B_o;
% T_rho{2,3} = 0;
% T_rho{3,1} = 0;
% T_rho{3,2} = 0;
% T_rho{3,3} = 1;


if any(computeLocalGas)
    q_g(computeLocalGas) =  B_g(computeLocalGas)./D(computeLocalGas).*      q_g_sc(computeLocalGas) - B_g(computeLocalGas).*R_s(computeLocalGas)./D(computeLocalGas).*   q_o_sc(computeLocalGas);
end
q_o = -B_o*r_s./D.*	q_g_sc + B_o./D.*       q_o_sc;
q_w =                                                  q_w_sc;

if any(computeLocalGas)
    rho_g(computeLocalGas) =  1./B_g(computeLocalGas).*	rho_sc(1) + r_s./B_g(computeLocalGas).*	rho_sc(2);
end
rho_o =  R_s./B_o.*	rho_sc(1) + 1./B_o.*	rho_sc(2);
rho_w =                                                rho_sc(3);




