function [dpds_tot,f_work,Z] = Beggs_Brill_dpds(s,p,notUsed,alpha,d,e,oil,q_g_sc,q_o_sc,q_w_sc,rho_sc,s_in,s_out,T_in,T_out,hasSurfaceGas,f_work,Z)
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
%Vectorization
%}


if nargin < 16
    hasSurfaceGas = true;
end
if nargin < 17 || isempty(f_work) || any(isnan(f_work))
    f_work = nan(size(double(p)));
end
if nargin < 18 || isempty(Z) || any(isnan(Z))
    Z = nan(size(double(p)));
end

if min(double(p)) < 1e5
    warning('Pressure below atmospheric.')
    dpds_tot  = 0; 
    return
end

zeroVal = (q_g_sc*0+q_o_sc*0+q_w_sc*0+p*0);  %% any of these can be an ADI

% % Determine inclination in case of survey file input:
% if length(alpha) > 1
%     n_sur = length(alpha(:,1)); % number of survey points
%     if s < alpha(1,1)
%         help = alpha(1,2);
%     else if s > alpha(n_sur,1)
%         help = alpha(n_sur,2);
%         else
%             help = interp1(alpha(:,1),alpha(:,2),s);
%         end
%     end
%     clear alpha;
%     alpha = help; % replace survey file by single inclination value, rad
% end

% Compute internal variables:
epsilon = e./d; % dimensionless pipe rougness, -
g = 9.81; % acceleration of gravity, m/s^2

% Compute local gas and liquid properties:
[mu_g,mu_l,q_g,q_l,rho_g,rho_l,sigma_gl,v_sg,v_sl,Z] = local_gas_liq_props(d,oil,p,...
    q_g_sc,q_o_sc,q_w_sc,rho_sc,s,s_in,s_out,T_in,T_out,hasSurfaceGas,Z);

% Check for free gas:
dpds_tot = zeroVal;
flow_reg = zeroVal;
cond_no_Local_gas_Flow =  abs(q_g) < 1.e-12;
if any(cond_no_Local_gas_Flow) % no free gas - liquid flow only
    
    [dpds_tot(cond_no_Local_gas_Flow),flow_reg(cond_no_Local_gas_Flow),f_work(cond_no_Local_gas_Flow)] = ...
        no_local_gas_flow(...
        v_sl(cond_no_Local_gas_Flow),...
        rho_l(cond_no_Local_gas_Flow),...
        d(cond_no_Local_gas_Flow),...
        mu_l(cond_no_Local_gas_Flow),...
        f_work(cond_no_Local_gas_Flow),...
        epsilon(cond_no_Local_gas_Flow),...
        g,...
        alpha(cond_no_Local_gas_Flow));

end
if any(~cond_no_Local_gas_Flow);
    [dpds_tot(~cond_no_Local_gas_Flow),flow_reg(~cond_no_Local_gas_Flow),f_work(~cond_no_Local_gas_Flow)] = ...
        NOT_cond_no_Local_gas_Flow(v_sl(~cond_no_Local_gas_Flow),...
        rho_l(~cond_no_Local_gas_Flow),...
        g,...
        sigma_gl(~cond_no_Local_gas_Flow),...
        v_sg(~cond_no_Local_gas_Flow),...
        q_l(~cond_no_Local_gas_Flow),...
        q_g(~cond_no_Local_gas_Flow),...
        alpha(~cond_no_Local_gas_Flow),...
        f_work(~cond_no_Local_gas_Flow),...
        d(~cond_no_Local_gas_Flow),...
        zeroVal(~cond_no_Local_gas_Flow),...
        rho_g(~cond_no_Local_gas_Flow),...
        mu_g(~cond_no_Local_gas_Flow),...
        mu_l(~cond_no_Local_gas_Flow),...
        epsilon(~cond_no_Local_gas_Flow),...
        p(~cond_no_Local_gas_Flow));
end

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

ee = ones(size(flow_dir));
ff = ones(size(flow_dir));
gg = ones(size(flow_dir));
hh = ones(size(flow_dir));

ee(flow_dir == 1 & flow_reg == 1) =  0.011;
ff(flow_dir == 1 & flow_reg == 1) = -3.7680;
gg(flow_dir == 1 & flow_reg == 1) =  3.5390;
hh(flow_dir == 1 & flow_reg == 1) = -1.6140;

ee(flow_dir == 1 & flow_reg == 2) =  2.960;
ff(flow_dir == 1 & flow_reg == 2) =  0.3050;
gg(flow_dir == 1 & flow_reg == 2) = -0.4473;
hh(flow_dir == 1 & flow_reg == 2) =  0.0978;

ee(flow_dir ~= 1) =  4.700;
ff(flow_dir ~= 1) = -0.3692;
gg(flow_dir ~= 1) =  0.1244;
hh(flow_dir ~= 1) = -0.5056;



C = (1.0-lambda_l).*log(ee.*lambda_l.^ff.*N_lv.^gg.*N_Fr.^hh);

C(C < 0) = 0;

C(flow_reg == 3 & flow_dir == 1) = 0;

Psi = 1.0 + C.*(sin(1.8*theta_BB)-0.333*(sin(1.8*theta_BB)).^3); % correction factor
H_l = H_l_0 .* Psi;
end


function [dpds_tot,flow_reg,f_work] = no_local_gas_flow(v_sl,rho_l,d,mu_l,f_work,epsilon,g,alpha)
    
    
    flow_reg = zeros(size(double(v_sl))); % liquid-only flow
    
    % Compute pressure gradient for liquid-only flow:
    v_l = v_sl; % local liquid velocity, m/s
    N_Re = rho_l.*d.*abs(v_l)./mu_l; % Reynolds number, -
    if isempty(f_work) || any(isnan(f_work))
        [f,f_work] = Moody_friction_factor(epsilon,N_Re); % friction factor, -
    else
        [f,f_work] = Moody_friction_factor(epsilon,N_Re,f_work); % friction factor, -
    end
    dpds_grav = rho_l.*g.*cos(alpha); % gravity losses, Pa/m
    dpds_fric = -rho_l.*f.*v_l.*abs(v_l)./(2.*d); % friction losses, Pa/m
    dpds_acc = 0; % acceleration losses are neglegible, Pa/m
    dpds_tot = dpds_grav + dpds_fric + dpds_acc; % total pressure gradient, Pa/m
    
end

function [dpds_tot,flow_reg,f_work] = NOT_cond_no_Local_gas_Flow(v_sl,rho_l,g,sigma_gl,v_sg,q_l,q_g,alpha,f_work,d,zeroVal,rho_g,mu_g,mu_l,epsilon,p)


    % Compute auxiliary variables:
    N_lv = abs(v_sl).*(rho_l./(g.*sigma_gl)).^(1/4); % liquid velocity number, -
    v_ms = v_sg + v_sl; % local mixture velocity, m/s
    lambda_l = q_l./(q_l+q_g); % local liquid volume fraction, -
    lambda_g = 1-lambda_l; % local gas volume fraction, -
    
    % Determine flow direction (uphill, downhill or horizontal)
    flow_dir = calculate_flow_dir(v_ms,alpha);

    
    % Determine the value of theta_BB. This is the angle as defined in the original publication
    % of Beggs and Brill.
    theta_BB = flow_dir.*abs(alpha-pi/2); % theta_BB is negative for downward and positive for
                                         % upward flow, rad  
  
    % Determine flow pattern as if pipe were horizontal:
    N_Fr = abs(v_ms).^2./(g*d);
    L_1 = 316*lambda_l.^0.302;
    L_2 = 0.000925*lambda_l.^-2.468;
    L_3 = 0.10*lambda_l.^-1.452;
    L_4 = 0.5*lambda_l.^-6.738;
    
    flow_reg = calculate_flow_reg(lambda_l,N_Fr,L_1,L_2,L_3,L_4);
   
     
    % Determine liquid hold-up as if pipe were horizontal:
    [aa,bb,cc] = getHoldupCoeficients(flow_reg);

    transitionFlow = flow_reg == 4;
    H_l_0 = zeroVal;
    if any(~transitionFlow) % outside transition flow regime
        H_l_0(~transitionFlow) = (aa(~transitionFlow).*lambda_l(~transitionFlow).^bb(~transitionFlow))./N_Fr(~transitionFlow).^cc(~transitionFlow);
    end
    if any(transitionFlow)% transition flow; interpolate
        aa(transitionFlow) =  0.980;
        bb(transitionFlow) =  0.4846;
        cc(transitionFlow) =  0.0868;
        H_l_0_seg = (aa(transitionFlow).*lambda_l(transitionFlow).^bb(transitionFlow))./N_Fr(transitionFlow).^cc(transitionFlow);
        aa(transitionFlow) =  0.845;
        bb(transitionFlow) =  0.5351;
        cc(transitionFlow) =  0.0173;
        H_l_0_int = (aa(transitionFlow).*lambda_l(transitionFlow).^bb(transitionFlow))./N_Fr(transitionFlow).^cc(transitionFlow);
        A = (L_3(transitionFlow)-N_Fr(transitionFlow))./(L_3(transitionFlow)-L_2(transitionFlow));
        H_l_0(transitionFlow) = A .* H_l_0_seg + (1-A) .* H_l_0_int;  
    end
    reality_check = H_l_0 < lambda_l;
    if any(reality_check) % reality check
        H_l_0(reality_check) = lambda_l(reality_check);
    end
    
    % Determine liquid hold up, corrected for pipe inclination:
    % Note: use is made of the local function Beggs_Brill_holdup, defined at the bottom of this
    % file.
    horizontalPipe = theta_BB == 0;
    H_l = H_l_0;  % set for all, then correct for those pipelines which are not horizontal
        
    if any(~horizontalPipe) % non-horizontal pipe        
        distuphill = flow_reg == 3 & flow_dir == 1;
        %if ~horizontalPipe & distuphill % distributed uphill flow; no correction
            %H_l(~horizontalPipe & distuphill) = H_l_0(~horizontalPipe & distuphill);
        if any(~horizontalPipe & ~distuphill)
            outsideTransition = flow_reg ~= 4;
            
            c1 = ~horizontalPipe & ~distuphill & outsideTransition;
            if any(c1) % outside transition flow regime   
                H_l(c1) = Beggs_Brill_holdup(flow_dir(c1),flow_reg(c1),H_l_0(c1),lambda_l(c1),N_Fr(c1),N_lv(c1),theta_BB(c1));
            end
            c2 =  ~horizontalPipe & ~distuphill & ~outsideTransition;
            if any(c2) % in transition flow regime; interpolate
                flow_reg(c2) = 1; 
                H_l_seg = Beggs_Brill_holdup(flow_dir(c2),flow_reg(c2),H_l_0(c2),lambda_l(c2),N_Fr(c2),N_lv(c2),...
                    theta_BB(c2));
                flow_reg(c2) = 2;
                H_l_int = Beggs_Brill_holdup(flow_dir(c2),flow_reg(c2),H_l_0(c2),lambda_l(c2),N_Fr(c2),N_lv(c2),...
                    theta_BB(c2));
                A = (L_3(c2)-N_Fr(c2))./(L_3(c2)-L_2(c2));
                H_l(c2) = A .* H_l_seg + (1-A) .* H_l_int;  
            end
        end
    end
    
    % Determine friction factor:
    rho_mn = rho_g.*lambda_g + rho_l.*lambda_l; % 'no-slip' mixture density, kg/m^3 
    mu_mn = mu_g.*lambda_g + mu_l.*lambda_l; % 'no-slip' mixture viscosity, Pa s
    rho_ms = rho_g.*(1-H_l) + rho_l.*H_l; % 'slip' mixture density, kg/m^3
    N_Re = rho_mn.*abs(v_ms).*d./mu_mn; % mixture Reynolds number, -
    y = lambda_l./H_l.^2;
    
    ss = y;  % initialization of the variable
    
    cy = (1 <= y & y <= 1.2);
    if any(cy)
        ss(cy) = log(2.2*y(cy)-1.2);
    end
    if any(~cy)
        ss(~cy) = log(y(~cy))./(-0.0523+3.182*log(y(~cy))-0.8725.*(log(y(~cy))).^2+0.01853.*(log(y(~cy))).^4);
    end
    
    f_f_n = exp(ss); % correction factor, -
    if isempty(f_work) || any(isnan(f_work))
        [f_n,f_work] = Moody_friction_factor(epsilon,N_Re); % 'normalising' friction factor, -
    else
        [f_n,f_work] = Moody_friction_factor(epsilon,N_Re,f_work); % 'normalising' friction factor, -
    end
    f = f_n.*(f_f_n); % friction factor, -
    
    % Determine pressure gradient:
    help01 = rho_mn.*v_ms.*v_sg./p; % acceleration loss factor E_k, -
    if any(help01 >= 1) 
        error('Choked flow. Reduce rate, increase diameter, or increase back pressure.')
    end
    help02 = rho_ms.*g.*cos(alpha);
    help03 = -f.*rho_mn.*v_ms.*abs(v_ms)./(2*d); 
    %dpds_grav = help02; % gravity losses, Pa/m
    %dpds_fric = help03; % friction losses, Pa/m
    dpds_tot = (help02+help03)./(1-help01); 
    %dpds_acc = help01*dpds_tot; % acceleration losses, Pa/m
    


end


function flow_dir = calculate_flow_dir(v_ms,alpha)

positive_vms = v_ms > 0;
downhill_well = alpha < pi/2;
uphill_well = alpha > pi/2;

flow_dir = nan(size(positive_vms));

flow_dir(positive_vms & downhill_well) = -1;
flow_dir(positive_vms & uphill_well) = 1;
flow_dir(positive_vms & ~downhill_well & ~uphill_well) = 0; %% horizontal

flow_dir(~positive_vms & downhill_well) = 1;
flow_dir(~positive_vms & uphill_well) = -1;
flow_dir(~positive_vms & ~downhill_well & ~uphill_well) = 0; %% horizontal

assert(~any(isnan(flow_dir)),'check flow direction classification')


end


function flow_reg = calculate_flow_reg(lambda_l,N_Fr,L_1,L_2,L_3,L_4)
flow_reg = nan(size(double(lambda_l)));

flow_reg1 = ((lambda_l < 0.01) & (N_Fr < L_1)) | ((lambda_l >= 0.01) & (N_Fr < L_2));

if any(flow_reg1)
    flow_reg(flow_reg1) = 1; % segregated flow
end
if any(~flow_reg1)
    flow_reg4 = (lambda_l >= 0.01 & L_2 <= N_Fr & N_Fr <= L_3);
    if any(~flow_reg1 & flow_reg4)
        flow_reg(~flow_reg1 & flow_reg4) = 4; % transition flow
    end
    if any(~flow_reg1 & ~flow_reg4)
        flow_reg2 = ((0.01 <= lambda_l & lambda_l < 0.4 & L_3 < N_Fr & N_Fr <= L_1) ...
                | (lambda_l >= 0.4 & L_3 < N_Fr & N_Fr <= L_4));
        if any(~flow_reg1 & ~flow_reg4 & flow_reg2)
            flow_reg(~flow_reg1 & ~flow_reg4 & flow_reg2) = 2;
        end
        if any(~flow_reg1 & ~flow_reg4 & ~flow_reg2)
            flow_reg(~flow_reg1 & ~flow_reg4 & ~flow_reg2) = 3;
        end
    
    end
    
end
assert(~any(isnan(flow_reg)),'check flow regime classification')

end

function [aa,bb,cc] = getHoldupCoeficients(flow_reg)

aa = zeros(size(flow_reg));
bb = zeros(size(flow_reg));
cc = zeros(size(flow_reg));

aa(flow_reg == 1) =  0.980;
bb(flow_reg == 1) =  0.4846;
cc(flow_reg == 1) =  0.0868;

aa(flow_reg == 2) =  0.845;
bb(flow_reg == 2) =  0.5351;
cc(flow_reg == 2) =  0.0173;

aa(flow_reg == 3) = 1.065;
bb(flow_reg == 3) = 0.5824;
cc(flow_reg == 3) = 0.0609;

end