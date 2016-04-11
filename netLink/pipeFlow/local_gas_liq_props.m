function [mu_g,mu_l,q_g,q_l,rho_g,rho_l,sigma_gl,v_sg,v_sl,Z] = local_gas_liq_props(d,oil,...
    p,q_g_sc,q_o_sc,q_w_sc,rho_sc,s,s_in,s_out,T_in,T_out,hasSurfaceGas,Z)
% function [mu_g,mu_l,q_g,q_l,rho_g,rho_l,sigma_gl,v_sg,v_sl] = local_gas_liq_props(d,oil,...
%     p,q_sc,rho_sc,s,s_in,s_out,T_in,T_out)
%
% Computes local gas and liquid properties for pipe flow.   
%
% alpha = inclination w.r.t. vertical, rad;
% d = inside diameter, m
% e = roughness, m
% p = pressure, Pa
% mu_g = local gas viscosity, Pa s
% mu_l = local liquid viscosity, Pa s
% oil = parameter to select black oil model or volatile oil table, -  
%   oil = 1: black oil; parameters computed with the aid of Standing correlations
%   oil = 2: black oil; parameters computed with the aid of Glaso correlations
%   oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
% q_g = local gas flow rate, m^3/s
% q_l = local liquid flow rate, m^3/s
% q_sc = [q_g_sc,q_o_sc,q_w_sc], m^3/s
%   q_g_sc = gas flow rate at standard conditions, m^3/s
%   q_o_sc = oil flow rate at standard conditions, m^3/s
%   q_w_sc = water flow rate at standard conditions, m^3/s
% rho_g = local gas density, kg/m^3
% rho_l = local liquid density, kg/m^3 
% rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc], kg/m^3 
%   rho_g_sc = gas density at standard conditions, kg/m^3
%   rho_o_sc = oil density at standard conditions, kg/m^3
%   rho_w_sc = water density at standard conditions, kg/m^3
% s = along-hole distance, measured from the separator to the reservoir, m   
% sigma_gl = local gas-liquid interfacial tension, N/m 
% s_in = starting point for the integration
% s_out = end point for the integration
% T_in = temperature at s_in, deg. C 
% T_out = temperature at s_out, deg. C
% v_sg = superficial gas velocity, m/s
% v_sl = superficial liquid velocity, m/s
%
% JDJ, 30-10-11, last revised 09-05-13
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%Vectorization
%Prevent unnecessary calculations when q_g_sc == 0
%}

if nargin < 13
    hasSurfaceGas = true;
end
if nargin < 14
    Z = [];
end

% Compute internal variables:
A = (pi*d.^2)./4; % pipe cross-sectional area, m^2

% Compute temperature through linear interpolation between T_in and T_out:
T = T_in+(T_out-T_in)*(s-s_in)./(s_out-s_in); % temperature, deg. C

% Densities and flow rates at standard conditions:
rho_g_sc = rho_sc(1); % gas density at standard conditions, kg/m^3
rho_o_sc = rho_sc(2); % oil density at standard conditions, kg/m^3
%q_g_sc = ; % gas flow rate at standard conditions, m^3/s
if ~hasSurfaceGas
    assert(all(double(q_g_sc)==0));
    q_g_sc = zeros(size(double(q_g_sc)));
end
%q_o_sc = q_sc(2); % oil flow rate at standard conditions, m^3/s

% Compute local gas and liquid properties:
if hasSurfaceGas
    R_go = q_g_sc./q_o_sc; % producing GOR as would be observed at surface, m^3/m^3
    R_sb = R_go; % This is the bubble point GOR for the oil in the wellbore. This value may be much
             % higher than R_sb in the reservoir if gas-cap gas or lift gas is produced.
else
    R_sb = zeros(size(q_g_sc));
end
[q_g,q_o,q_w,rho_g,rho_o,rho_w,Z] = local_q_and_rho(oil,p,q_g_sc,q_o_sc,q_w_sc,R_sb,rho_sc,T,hasSurfaceGas,Z);
% q = [q_g, q_o, q_w], m^3/s, rho = [rho_g, rho_o, rho_w], kg/m^3 

%q_g = q(1); % local gas flow rate, m^3/s
%q_o = q(2); % local oil flow rate, m^3/s
%q_w = q(3); % local water flow rate, m^3/s

%rho_g = rho(1); % local gas density, kg/m^3
%rho_o = rho(2); % local oil density, kg/m^3
%rho_w = rho(3); % local water density, kg/m^3

c1 = abs(double(q_g)) >= 1e-12;
mu_g = p;
if any(c1)
    mu_g(c1) = gas_viscosity(p(c1),rho_g_sc,T(c1)); % local gas viscosity, Pa s
end
if any(~c1)
    mu_g(~c1) = 0;    
end
mu_o = oil_viscosity(p,R_sb,rho_g_sc,rho_o_sc,T); % local oil viscosity, Pa s
mu_w = water_viscosity; % input function; local water viscosity, Pa s

sigma = interfacial_tensions; % input function; sigma = [sigma_go, sigma_gw]; 
sigma_go = sigma(1); % gas-oil interfacial tension, N/m
sigma_gw = sigma(2); % gas-water interfacial tension, N/m

f_o = q_o./(q_o+q_w); % local oil fraction , - 
f_w = q_w./(q_o+q_w); % local water fraction, -

q_l = q_o + q_w; % local liquid flow rate, m^3/s 
rho_l = rho_o.*f_o + rho_w.*f_w; % local liquid density, kg/m^3 
mu_l = mu_o.*f_o + mu_w.*f_w; % local liquid viscosity, Pa s
sigma_gl = sigma_go*f_o + sigma_gw*f_w; % local gas-liquid interfacial tension, N/m 

% Compute superficial velocities:
v_sg = q_g./A; % local superficial gas velocity, m/s 
v_sl = q_l./A; % local superficial liquid velocity, m/s


