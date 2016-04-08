function dpds = Beggs_Brill_dpds(s,p,~,alpha,d,e,oil,q_sc,rho_sc,s_in,s_out,T_in,T_out)
% function dpds = Beggs_Brill_dpds(s,p,~,alpha,d,e,oil,q_sc,rho_sc,s_in,s_out,T_in,T_out)
%
% Computes the derivative dp/ds for a given pressure p and along-hole distance s, in an element
% of a flowline-wellbore system. The distance s is measured from the separator towards the
% reservoir. Therefore, flowrates are negative for production wells.
%
% Uses the Beggs and Brill correlation for multiphase flow in inclined wells; see references
% [1] and [2]. Corrections have been implemented as suggested by Payne et al. [2], [3].        
%
% The vector p contains the total pressure, and the pressures taking into account the
% individual effects of gravity, friction and acceleration losses respectively. Accordingly,
% the vector dpds contains the total pressure loss per unit length, as well as the individual
% gravity losses, friction losses and acceleration losses.
%
% This function can be used to compute the pressure drop through numerical integration. It has
% the correct format to be used in conjunction with one of the standard numerical integration
% routines in MATLAB.    
%
% alpha = inclination wrt. vertical, rad; alternatively alpha can be a survey file (matrix)
%         with AHD values in the first column (in m) and inclination values in the second
%         column (in rad).    
% d = inside diameter, m
% dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc]
%   dpds_acc =  pressure gradient due to acceleration losses, Pa/m
%   dpds_fric =  pressure gradient due to friction losses, Pa/m 
%   dpds_grav = pressure gradient due to head losses, Pa/m 
%   dpds_tot = dpds_grav + dpds_fric + dpds_acc = total pressure gradient, Pa/m  
% e = roughness, m
% oil = parameter to select black oil model or volatile oil table, -  
%   oil = 1: black oil; parameters computed with the aid of Standing correlations
%   oil = 2: black oil; parameters computed with the aid of Glaso correlations
%   oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
% p = [p_tot,p_grav,p_fric,p_acc], Pa
%   p_acc  = p_in + pressure increase (decrease for production wells) due to accel. losses, Pa
%   p_fric = p_in + pressure increase (decrease for production wells) due to frict. losses, Pa
%   p_grav = p_in + pressure increase (decrease for production wells) due to head loss, Pa
%   p_tot  = p_in + pressure increase (decrease for production wells) due to gravity, friction
%            and acceleration losses, Pa 
% p_in = pressure at s_in, Pa 
% p_out = pressure at s_out, Pa
% q_sc = [q_g_sc,q_o_sc,q_w_sc], m^3/s
%   q_g_sc = gas flow rate at standard conditions, m^3/s
%   q_o_sc = oil flow rate at standard conditions, m^3/s
%   q_w_sc = water flow rate at standard conditions, m^3/s
% rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc], kg/m^3 
%   rho_g_sc = gas density at standard conditions, kg/m^3
%   rho_o_sc = oil density at standard conditions, kg/m^3
%   rho_w_sc = water density at standard conditions, kg/m^3
% s = along-hole distance, measured from the separator to the reservoir, m   
% s_in = starting point for the integration
% s_out = end point for the integration
% T_in = temperature at s_in, deg. C 
% T_out = temperature at s_out, deg. C 
%
% JDJ, JHvB, 26-09-03, last revised JDJ 10-05-13
%
% References:
% [1] Beggs, H.D. and Brill, J.P., 1973: A study of two-phase flow in inclined pipes, J. Petr.
%     Techn., May, p.607-617.  
% [2] Brill, J.P. and Mukherjee, H., 1999: Multiphase flow in wells, SPE Monograph Series,
%     vol. 17, SPE, Richardson. 
% [3] Payne, G.A., Palmer, C.M., Brill, J.P. and Beggs, H.D., 1979: Evaluation of
%     inclined-pipe, two-phase liquid holdup and pressure-loss correlations using experimental
%     data, J. Petr. Techn., September, p. 1198-1208.     

% Check sign of pressure:
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%}






p_tot = p(1); % first element of vector p is the total wellbore pressure, Pa 
if min(double(p_tot)) < 1e5
    warning('Pressure below atmospheric.')
    dpds_tot  = 0; 
    dpds_grav = 0; 
    dpds_fric = 0; 
    dpds_acc  = 0; 
    dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc];
    return
end

% Determine inclination in case of survey file input:
if length(alpha) > 1
    n_sur = length(alpha(:,1)); % number of survey points
    if s < alpha(1,1)
        help = alpha(1,2);
    else if s > alpha(n_sur,1)
        help = alpha(n_sur,2);
        else
            help = interp1(alpha(:,1),alpha(:,2),s);
        end
    end
    clear alpha;
    alpha = help; % replace survey file by single inclination value, rad
end

% Compute internal variables:
epsilon = e/d; % dimensionless pipe rougness, -
g = 9.81; % acceleration of gravity, m/s^2

% Compute local gas and liquid properties:
[mu_g,mu_l,q_g,q_l,rho_g,rho_l,sigma_gl,v_sg,v_sl] = local_gas_liq_props(d,oil,p_tot,...
    q_sc,rho_sc,s,s_in,s_out,T_in,T_out);

% Check for free gas:
if abs(q_g) < 1.e-12 % no free gas - liquid flow only
    flow_reg = 0; % liquid-only flow
    
    % Compute pressure gradient for liquid-only flow:
    v_l = v_sl; % local liquid velocity, m/s
    N_Re = rho_l*d*abs(v_l)./mu_l; % Reynolds number, -
    f = Moody_friction_factor(epsilon,N_Re); % friction factor, -
    dpds_grav = rho_l*g*cos(alpha); % gravity losses, Pa/m
    dpds_fric = -rho_l*f*v_l*abs(v_l)./(2*d); % friction losses, Pa/m
    dpds_acc = 0; % acceleration losses are neglegible, Pa/m
    dpds_tot = dpds_grav + dpds_fric + dpds_acc; % total pressure gradient, Pa/m
    dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc];
    
else % gas-liquid flow

    % Compute auxiliary variables:
    N_lv = abs(v_sl)*(rho_l./(g*sigma_gl)).^(1/4); % liquid velocity number, -
    v_ms = v_sg + v_sl; % local mixture velocity, m/s
    lambda_l = q_l./(q_l+q_g); % local liquid volume fraction, -
    lambda_g = 1-lambda_l; % local gas volume fraction, -
    
    % Determine flow direction (uphill, downhill or horizontal)
    if v_ms > 0 % flow from wellhead to bottomhole (injection well)
        if  alpha < pi/2 % 'downhill' drilled well section (usual situation)
            flow_dir = -1; % downhill flow
        else
            if  alpha > pi/2 % 'uphill' drilled well section (occurs
                             % occasionally in 'horizontal' wells)
                flow_dir = 1; % uphill flow
            else % alpha = pi/2, horizontal well section
                flow_dir = 0; % horizontal flow
            end
        end
    else % flow from bottomhole to wellhead (production well)
        if  alpha < pi/2 % 'downhill' drilled well section
            flow_dir = 1; % uphill flow
        else
            if alpha > pi/2 % 'uphill' drilled well section
                flow_dir = -1; % downhill flow
            else % alpha = pi/2, horizontal well section
                flow_dir = 0; % horizontal flow
            end
        end
    end
    
    % Determine the value of theta_BB. This is the angle as defined in the original publication
    % of Beggs and Brill.
    theta_BB = flow_dir*abs(alpha-pi/2); % theta_BB is negative for downward and positive for
                                         % upward flow, rad  
  
    % Determine flow pattern as if pipe were horizontal:
    N_Fr = abs(v_ms).^2./(g*d);
    L_1 = 316*lambda_l.^0.302;
    L_2 = 0.000925*lambda_l.^-2.468;
    L_3 = 0.10*lambda_l.^-1.452;
    L_4 = 0.5*lambda_l.^-6.738;
    if ((lambda_l < 0.01) && (N_Fr < L_1)) || ((lambda_l >= 0.01) && (N_Fr < L_2)) 
        flow_reg = 1;  % segregated flow
    else
        if (lambda_l >= 0.01 && L_2 <= N_Fr && N_Fr <= L_3) 
            flow_reg = 4; % transition flow
        else
            if ((0.01 <= lambda_l && lambda_l < 0.4 && L_3 < N_Fr && N_Fr <= L_1) ...
                    || (lambda_l >= 0.4 && L_3 < N_Fr && N_Fr <= L_4))
                flow_reg = 2; % intermittent flow
            else 
                flow_reg = 3; % distributed flow
            end
        end
    end    
     
    % Determine liquid hold-up as if pipe were horizontal:
    if flow_reg == 1 % segregated flow
        aa =  0.980;
        bb =  0.4846;
        cc =  0.0868;
    else
        if flow_reg == 2 % intermittent flow
            aa =  0.845;
            bb =  0.5351;
            cc =  0.0173;
        else
            if flow_reg == 3 % distributed flow
                aa = 1.065;
                bb = 0.5824;
                cc = 0.0609;              
            end
        end
    end
    if flow_reg ~= 4 % outside transition flow regime
        H_l_0 = (aa*lambda_l.^bb)./N_Fr.^cc;
    else % transition flow; interpolate
        aa =  0.980;
        bb =  0.4846;
        cc =  0.0868;
        H_l_0_seg = (aa*lambda_l.^bb)./N_Fr.^cc;
        aa =  0.845;
        bb =  0.5351;
        cc =  0.0173;
        H_l_0_int = (aa*lambda_l.^bb)./N_Fr.^cc;
        A = (L_3-N_Fr)./(L_3-L_2);
        H_l_0 = A * H_l_0_seg + (1-A) * H_l_0_int;  
    end
    if H_l_0 < lambda_l % reality check
        H_l_0 = lambda_l;
    end
    
    % Determine liquid hold up, corrected for pipe inclination:
    % Note: use is made of the local function Beggs_Brill_holdup, defined at the bottom of this
    % file.
    if theta_BB == 0 % horizontal pipe
        H_l = H_l_0;
    else % non-horizontal pipe        
        if flow_reg == 3 && flow_dir == 1 % distributed uphill flow; no correction
            H_l = H_l_0;
        else
            if flow_reg ~= 4 % outside transition flow regime   
                H_l = Beggs_Brill_holdup(flow_dir,flow_reg,H_l_0,lambda_l,N_Fr,N_lv,theta_BB);
            else % in transition flow regime; interpolate
                flow_reg = 1; 
                H_l_seg = Beggs_Brill_holdup(flow_dir,flow_reg,H_l_0,lambda_l,N_Fr,N_lv,...
                    theta_BB);
                flow_reg = 2;
                H_l_int = Beggs_Brill_holdup(flow_dir,flow_reg,H_l_0,lambda_l,N_Fr,N_lv,...
                    theta_BB);
                A = (L_3-N_Fr)./(L_3-L_2);
                H_l = A * H_l_seg + (1-A) * H_l_int;  
            end
        end
    end
    
    % Determine friction factor:
    rho_mn = rho_g*lambda_g + rho_l*lambda_l; % 'no-slip' mixture density, kg/m^3 
    mu_mn = mu_g*lambda_g + mu_l*lambda_l; % 'no-slip' mixture viscosity, Pa s
    rho_ms = rho_g*(1-H_l) + rho_l*H_l; % 'slip' mixture density, kg/m^3
    N_Re = rho_mn*abs(v_ms)*d./mu_mn; % mixture Reynolds number, -
    y = lambda_l./H_l.^2;
    if (1 <= y && y <= 1.2)
        ss = log(2.2*y-1.2);
    else
        ss = log(y)./(-0.0523+3.182*log(y)-0.8725*(log(y)).^2+0.01853*(log(y)).^4);
    end
    f_f_n = exp(ss); % correction factor, -
    f_n = Moody_friction_factor(epsilon,N_Re); % 'normalising' friction factor, -
    f = f_n*(f_f_n); % friction factor, -
    
    % Determine pressure gradient:
    help01 = rho_mn*v_ms*v_sg./p_tot; % acceleration loss factor E_k, -
    if help01 >= 1 
        error('Choked flow. Reduce rate, increase diameter, or increase back pressure.')
    end
    help02 = rho_ms*g*cos(alpha);
    help03 = -f*rho_mn*v_ms*abs(v_ms)./(2*d); 
    dpds_grav = help02; % gravity losses, Pa/m
    dpds_fric = help03; % friction losses, Pa/m
    dpds_tot = (help02+help03)./(1-help01); 
    dpds_acc = help01*dpds_tot; % acceleration losses, Pa/m
	dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc];
    
end  

function H_l = Beggs_Brill_holdup(flow_dir,flow_reg,H_l_0,lambda_l,N_Fr,N_lv,theta_BB)
% function H_l = Beggs_Brill_holdup(flow_dir,flow_reg,H_l_0,lambda_l,N_Fr,N_lv,theta_BB)
%
% Determines liquid hold up, corrected for pipe inclination. Local function, used only in
% Beggs_Brill_dpds. 
%
% flow_dir = flag for flow direction, -
% flow_ regime = flag for flow regime, -
% H_l_0 = liquid holdup as if pipe were horizontal, -
% lambda_l = liquid volume fraction, -
% N_Fr = Froude number, -  
% N_lv = liquid velocity number, -
% theta_BB = inclination as defined in Beggs and Brill (1973), rad 
%
% JDJ, 26-09-03, last revised 10-05-13

if flow_dir == 1 % uphill flow
    if flow_reg == 1 % segregated flow
        ee =  0.011;
        ff = -3.7680;
        gg =  3.5390;
        hh = -1.6140;
    else
        if flow_reg == 2 % intermittent flow
            ee =  2.960;
            ff =  0.3050;
            gg = -0.4473;
            hh =  0.0978; 
        end
    end
else % downhill flow
    ee =  4.700;
    ff = -0.3692;
    gg =  0.1244;
    hh = -0.5056;
end
C = (1.0-lambda_l)*log(ee*lambda_l.^ff*N_lv.^gg*N_Fr.^hh);
if C < 0
    C = 0;
end
if flow_reg == 3 && flow_dir == 1 % distributed uphill flow
    C = 0;
end
Psi = 1.0 + C*(sin(1.8*theta_BB)-0.333*(sin(1.8*theta_BB)).^3); % correction factor
H_l = H_l_0 * Psi;    
