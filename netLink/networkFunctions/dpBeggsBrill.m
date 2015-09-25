% function dp = dpin_uphill(e, v)
%     %% calculates pressure drop in edge e from flow leaving vertex v.
% 
%     %% stream mock object
%     str = getStream(e);
%              
%     %% TODO: change this function to receive the inlet pressure
%     dp = dpBeggsBrill(e, str, v.pressure);
% end
% 
% function dp = dpout_downhill(e,v)
%     %% calculates pressure drop in edge e from flow arriving at v.
% 
%      %% stream mock object
%      str = getStream(e);
%      
%      dpBeggsBrill(e, str, v.pressure);
% end
% 
% function dp = dpin_downhill(e, v)
%     %% calculates pressure drop in edge e from flow leaving vertex v.
% 
%     %% stream mock object
%     str = getStream(e);
%              
%     %% TODO: change this function to receive the inlet pressure
%     dp = dpBeggsBrill(e, str, v.pressure);
% end



function dp_psi_tot =  dpBeggsBrill(E, qoE, qwE, qgE, pV)
%% dp: calculates the pressure drop for inlet flow in psi
%% TODO: extend this function to a vector of edges E = [e1, e2,.. e_m]
     
    units = struct('METRIC',0, 'FIELD', 1);    
    flow_regime = struct('SEGREGATED', 0, 'TRANSITION', 1, 'INTERMITTENT', 2, 'DISTRIBUTED', 3, 'UNDEFINED', 4);    
    
    str = vertcat(E.stream);       
    
    %% Producer flows are negative, and injection flows are positive.
    total_rate = qgE + qoE + qwE;   
    flag_rate = total_rate ~= 0;
    
    if (~flag_rate)
        dp_psi_tot = zeros(length(total_rate),1);
        return;
    end       
 
    %p_psi       = pV .* 14.5037788;                                % pressure in psi        
    p_psi = pV;
    
    
    pipes = vertcat(E.pipeline);
    diameters = vertcat(pipes.diam);
    temperatures = vertcat(pipes.temp);
    angles = vertcat(pipes.ang);
    
    diam_in     = diameters.*39.3700787;                       % pipe diameter in inches
    temp_f      = convtemp(temperatures,'C', 'F');              % temperature in F
    
    g           = 32.2;                                          % gravitational constant (ft/s^2)
    zfactor     = gasZFactor(vertcat(str.sg_gas), temperatures, p_psi);     % gas z-factor
    
    %%TODO: correct calculation of zfactor
    %zfactor = 0.935; (result = 0.97)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Flow Velocities %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vsl         = superficialLiquidVelocity(qoE, qwE,diam_in);                % superficial liquid velocity
    vsg         = superficialGasVelocity(qgE, p_psi, zfactor, diam_in, temp_f);                   % superficial two phase velocity
    vm          = vsl + vsg;                                     % superficial two phase velocity
    
    
    froude_num = vm.^2./diam_in./g;                                 % froude number    
    liquid_content = vsl./vm;
    
    
    %% if the liquid content is 0, changing it to something small
    liquid_content = max(liquid_content,  1e-8.*ones(length(E),1));
        
    %% constants used in formula for evaluating the flow regime  
    l1 = 316.*liquid_content.^0.302;
    l2 = 0.0009252.*liquid_content.^-2.4684;
    l3 = 0.10.*liquid_content.^(-1.4516);
    l4 = 0.5.*liquid_content.^(-6.738);
        
    
    %% checking the flow regime   
    %if (liquid_content < 0.01 & froude_num < l1) regime = flow_regime.SEGREGATED;
    %elseif(liquid_content >= 0.01 & froude_num < l2) regime = flow_regime.SEGREGATED;
    %elseif(liquid_content >= 0.01 & froude_num >= l2 & froude_no <= l3) regime = flow_regime.TRANSITION;
    %elseif(liquid_content >= 0.01 & liquid_content < 0.4 & froude_num > l3 & froude_num <= l1) regime = flow_regime.INTERMITTENT;
    %elseif(liquid_content >= 0.4 & froude_num > l3 & froude_num <= l4) regime = flow_regime.INTERMITTENT;
    %elseif(liquid_content < 0.4 & froude_num >= l1) regime = flow_regime.DISTRIBUTED;
    %elseif(liquid_content >= 0.4 & froude_num > l4) regime = flow_regime.DISTRIBUTED;
    %else
    %    disp '### Warning ###';
    %    disp 'From: Beggs & Brill 1973';
    %    disp 'The flow regime could not be determined';
    %    disp 'Assuming INTERMITTENT flow for the current stream';        
    %    
    %   	regime = flow_regime.INTERMITTENT;
    %end
    
    %% conditional vectors which are used to determine flow regimes
    cond_seg = (liquid_content < 0.01 & froude_num < l1) | (liquid_content>=0.01 & froude_num < l2);
    cond_trans = (liquid_content >= 0.01 & froude_num >= l2 & froude_num <= l3);    
    cond_dist = (liquid_content < 0.4 & froude_num >= l1) | (liquid_content >= 0.4 & froude_num > l4);
    cond_interm = (liquid_content >= 0.01 & liquid_content < 0.4 & froude_num > l3 & froude_num <= l1) | ...
                  (liquid_content >= 0.4 & froude_num > l3 & froude_num <= l4) | ...
                  (~cond_seg & ~cond_trans & ~cond_dist);  % assuming intermittent flow when it is not possible to determine.

              
    regime = zeros(length(E),1);
    
    %% Horizontal liquid holdup, hl(0), except for TRANSITION flow, 
    %% which will be later calculated as an interpolation of SEGREGATED and
    %% INTERMITTENT flows.
    hz_holdup = zeros(length(E),1);  
    
    if any(cond_seg) % Segregated flow
        regime(cond_seg) =  flow_regime.SEGREGATED; 
        hz_holdup(cond_seg) = (0.98.*(liquid_content(cond_seg).^0.4846))./(froude_num(cond_seg).^0.0868);  
    elseif any(cond_trans) % Transition flow
        regime(cond_trans) = flow_regime.TRANSITION;
    elseif any(cond_dist) % Distributed flow
        regime(cond_dist) = flow_regime.DISTRIBUTED;
        hz_holdup(cond_dist) = (1.065.*(liquid_content(cond_dist).^0.5824))./(froude_num(cond_dist).^0.0609);
    elseif any(cond_interm) % Intermittent flow    
        regime(cond_interm) = flow_regime.INTERMITTENT;
        hz_holdup(cond_interm) =  (0.845.*(liquid_content(cond_interm).^0.5351))./(froude_num(cond_interm).^0.0173);
    end
    
    assert(all(or(cond_interm,or(cond_seg,or(cond_trans,cond_dist)))),'Couldn''t find flow regime');        
   
   
    %% the horizontal holdup is equal to liquid content if smaller
    hz_holdup = max(hz_holdup, liquid_content);
   
   
   %% calculation of the correction factor   
   %TODO: verify next functions
   den_l = liquidDensity(qoE, qwE, str);                          %% liquid density
   den_g = gasDensity(gasSpecificGravity(str), temperatures, p_psi, zfactor);    %% gas density
   surface_tens = surfaceTension(den_g, den_l);          %% gas - liquid surface tension
   
   nlv = vsl .* (den_l./(g.*surface_tens)).^(0.25);         %% liquid velocity number
   
%    cor = 0.0;
%    payne_cor = 0.924;
%     if (angles < 0) 
%        payne_cor = 0.685;
%        cor = (1-liquid_content).*log(4.7.*(liquid_content.^-0.3692).*(nlv.^0.1244).*(froude_num.^-0.5056));
%    elseif (regime == flow_regime.SEGREGATED)
%        cor = (1 - liquid_content).*log(0.011.*(liquid_content.^-3.768).*(nlv.^3.539).*(froude_num.^-1.614));
%    elseif (regime == flow_regime.INTERMITTENT)
%        cor = (1 - liquid_content).*log(2.96 .*(liquid_content.^-0.305).*(nlv.^-0.4473).*(froude_num.^-0.0978));
%    elseif (regime == flow_regime.DISTRIBUTED)
%        cor = 0.0;
%    end

   
    %% Payne correction factor to holdup
   cor = zeros(length(E),1);                                 
   payne_cor = 0.924.*ones(length(E),1);   
  
   condAng = (angles < zeros(length(angles),1));
   if any(condAng) % Negative Angle       
       payne_cor(condAng) = 0.685;
       cor(condAng) = (ones(sum(condAng),1)-liquid_content(condAng)).*log(4.7.*(liquid_content(condAng).^-0.3692).*(nlv(condAng).^0.1244).*(froude_num(condAng).^-0.5056));
   elseif any(cond_seg) % SEGREGATED Flow       
       cor(cond_seg) = (ones(sum(cond_seg),1) - liquid_content(cond_seg)).*log(0.011.*(liquid_content(cond_seg).^-3.768).*(nlv(cond_seg).^3.539).*(froude_num(cond_seg).^-1.614));
   elseif any(cond_interm) % INTERMITTENT Flow    
       cor(cond_interm) = (ones(sum(cond_interm),1) - liquid_content(cond_interm)).*log(2.96 .*(liquid_content(cond_interm).^-0.305).*(nlv(cond_interm).^-0.4473).*(froude_num(cond_interm).^-0.0978));
   elseif any(cond_dist) % DISTRIBUTED flow       
       cor(cond_dist) = zeros(sum(cond_dist),1);
   end
   

   %% checking if correction >= 0
   condCor = (cor < 0);
   if any(condCor)
       disp '### Warning ###'
       disp 'From: Beggs & Brill 1973'
       disp 'The calculated correction factor, C, is negative... '
       disp 'Resetting to 0.0'
       
       cor(condCor) = zeros(sum(condCor),1);
   end
   
   %% liquid holdup correction for inclination
   phi = ones(length(angles),1) + cor.*(sin(pi.*1.8.*angles./180) - 0.333.*(sin(pi.*1.8.*angles./180).^3));
   holdup = payne_cor.*hz_holdup.*phi;
   
   
   %% if trasition regime, the liquid holdup is a mix of segregated and intermmitent
   if any(cond_trans)
       holdup(cond_trans) = liqholdupTransitionFlow(liquid_content(cond_trans), nlv(cond_trans), froude_num(cond_trans), angles(cond_trans), l2(cond_trans), l3(cond_trans));
   end
   
   %% corrects the holdup if wrong
   holdup = -(max(-holdup, -1.0.*ones(length(holdup),1)));
   
   
   %% calculating pressure drop due to elevation change
   den_s = den_l .* holdup + den_g .* (1 - holdup);  %% two phase density
   dp_el = den_s .* sin((pi.*angles)./180)./144;      %% pressure drop due to elevation change


   %% calculating friction factor    
   vis_g = gasViscosity(str, temperatures, p_psi, zfactor);
   den_ns = den_l .* liquid_content + den_g .* (1 - liquid_content);      %% no-slip density
   vis_ns = liquidViscosity(qoE, qwE, str) .* liquid_content + vis_g .* (1 - liquid_content);       %% no-slip viscosity   
   
   re_ns = 124.*(den_ns .* vm .* diam_in)./vis_ns;     %% no-slip reynolds number

   
   %%TODO: why denominator is log(-x) ? This is producing an irrational
   %%number.   
   fn = 1./(2 .* log((re_ns ./ (4.5223 .* log(re_ns)./log(10) - 3.8215)))./log(10)).^2;    %% no-slip friction factor
   
   %fn = 0.0056 + 0.5./(re_ns).^(0.32); % simple calculation for fn

   y = liquid_content ./(holdup.^2);    
   
   cond_y = (y > 1.0 & y < 1.2);
   s_term = log(y(~cond_y)) ./ (-0.0523 + 3.182 .* log(y(~cond_y)) - 0.8725 .* (log(y(~cond_y)).^2) + 0.01853 .* (log(y(~cond_y)).^4));
   if any(cond_y)
       s_term(cond_y) = log(2.2.*y(cond_y) -1.2); 
   end

   ftp = fn .* exp(s_term);      %% the friction factor

   %% calculating pressure drop due to friction   
   dp_f = 5.176e-3 .*(ftp .* den_ns .* (vm.^2)) ./ (diam_in);


   %% calculating acceleration term
   ek = 2.16e-4 .* (den_ns .* vm .* vsg) ./ p_psi;

   %% calculating total pressure drop
   
   dp_tot = (dp_f + dp_el) ./ (1 - ek);  %% total pressure drop per length of pipe (in psi ./ ft)
   %% double dp_tot_bar = dp_tot ./ 14.5038 ./ 0.3048;       // total pressure drop in bar / m


   %% converting length from m to ft
   length_ft =  vertcat(pipes.len).*3.28;  
   

   %% total pressure drop in psi
   dp_psi_tot = dp_tot .* length_ft;    
   
   cond_metric = (vertcat(E.units) == units.METRIC.*ones(length(E),1));
   if any(cond_metric)
       dp_psi_tot(cond_metric) = dp_psi_tot(cond_metric)./14.5037738;
   end 
end

function [holdup] = liqholdupTransitionFlow(liquid_content, nlv, froud_num, angles, l2, l3)
       frac = (l3-froud_num)./(l3-l2);
       
       %% horizontal holdups
       hz_holdup_seg = frac.*(0.98.*(liquid_content.^0.4846))./(froude_num.^0.0868);
       hz_holdup_int = (1-frac).*(0.845.*(liquid_content.^0.5351))./(froude_num.^0.0173);
       
       %% horizontal holdup is the liquid content if smaller.
       %% -(max(-x,-y)) = min(x,y) for positive integers
       hz_holdup_seg = -(max(-hz_holdup_seg, -liquid_content));       
       hz_holdup_int = -(max(-hz_holdup_int, -liquid_content));
       
       %% correction factors
       cor_seg = zeros(length(E),1);       
       cor_int = zeros(length(E),1);
       
       cond_ang = (angles < 0);
       if any(cond_ang)
           cor_seg(cond_ang) = (1-liquid_content(cond_ang)).*log(4.7.*(liquid_content(cond_ang).^-0.3692).*(nlv(cond_ang).^0.1244).*(froude_num(cond_ang).^-0.5056));
           cor_int(cond_ang) = cor_seg(cond_ang);
       else
           cor_seg(~cond_ang) = (1 - liquid_content(~cond_ang)) .* log(0.011 .*(liquid_content(~cond_ang).^-3.768).*(nlv(~cond_ang).^3.539).*(froude_num(~cond_ang).^-1.614));
           cor_int(~cond_ang) = (1 - liquid_content(~cond_ang)) .* log(2.96 .* (liquid_content(~cond_ang).^-0.305) .* (nlv(~cond_ang).^-0.4473) .*(froude_num(~cond_ang).^-0.0978));
       end 
       
       phi_seg = 1 + cor_seg .* (sin(pi .* 1.8.*angles./180) - 0.333 .* pow(sin(pi .* 1.8.*angles./ 180), 3));
       phi_int = 1 + cor_int .* (sin(pi .* 1.8.*angles./180) - 0.333 .* pow(sin(pi .* 1.8.*angles./ 180), 3));

       holdup = payne_cor .* (frac .* (hz_holdup_seg .* phi_seg) + (1 - frac) .* (hz_holdup_int .* phi_int));

end

%% gasSpecificGravity: returns the gas specific gravity
function sg = gasSpecificGravity(str)
    sg = vertcat(str.sg_gas);
end

% superficialGasVelocity: calculates the superficial gas velocity (vsg)
function gv = superficialGasVelocity(qgE, pres, zfac, diam, temp)    
    gas_rate_surface = qgE./0.0283168466;      % the gas rate in Sft^3/s
    
    % surface conditions
    Tsc = convtemp(15.65,'C','R');    
    Psc = 14.7; 
    
    % pipe conditions
    temp_R = convtemp(temp,'F','R');
                                                                                                                    
    A = pi.*((diam./2).*0.0833333).^2;     % pipe area                                                                        in ft
    gv = (gas_rate_surface.*zfac.*temp_R.*Psc)./(A.*Tsc.*pres.*86400);
end


%% superficialLiquidVelocity: calculates the superficial liquid velocity (vsl)
function liqRateFt =  superficialLiquidVelocity(qo, qw, diam)
    
    liquid_rate = 6.29.*(qo + qw);              % liquid rate in bbl/d
    liquid_rate_ft = 5.61458333.*liquid_rate./86400;            % liquid rate in ft^3/s     
                                
    
    den = pi.*((diam./2).*0.0833333).^2;     % pipe radius in ft
    
    liqRateFt = liquid_rate_ft./den;                 % in ft/s
end

%% liqDens: calculates the liquid density (oil and water)
function liqDens = liquidDensity(qoE, qwE, str)
    oil_rate = qoE;
    water_rate = qwE;
   
    if (oil_rate + water_rate) < 1.e-6
        oil_rate = ones(length(oil_rate),1);
    end
    
    
    oil_mdensity = oil_rate.*vertcat(str.oil_dens);
    water_mdensity =  water_rate.*vertcat(str.water_dens);
    
    
    den_metric = (oil_mdensity + water_mdensity)./(oil_rate + water_rate);   % liquid density
    liqDens =  0.0624279606.*den_metric;                                                     % converting to lb/ft^3        
end
    

function zfac = gasZFactor(sg, t, p)
%% gasZFactor: calculates gas z-factor
% t  -  oC
% p  - bara
% sg - gas specific gravity

%%TODO: evalute the partial differential df/dp (take into account z(p))
 
    %% calculating pseudo-critical properties (Shutton Equation)
    t_pc = 169.2 + 349.5.*sg - 74.*sg.^2; 
    p_pc = 756.8 - 131.*sg - 3.6.*sg.^2;    
 
    %% unit conversion
    t = convtemp(t,'C','F');                %% oC to F
    
    %% calculating pseud reduced properties
    t_pr = (t+460)./t_pc;
    p_pr = p./p_pc;   
    
    t = 1./t_pr;    
    a =  0.06125.* t .* (-1.2.*(1-t)).^2;
    
    y = 0.001.*ones(length(sg),1);
    i = 0;
    fy = ones(length(sg),1); % remove this line
    cond_fy =  abs(fy) > 1e-08;
    
    while true        
         [fy, dfy] = estimateZfactor(t(cond_fy), y(cond_fy), a, p_pr);                
         
         y(cond_fy) =  y(cond_fy) -(fy(cond_fy)./dfy(cond_fy)); 
        
        cond_fy = abs(fy) > 1e-08;
        
        if (all(~cond_fy) || i >= 200)
            break;     
        end
        
        i =i + 1;
    end    
    zfac = a.*p_pr./ y;
end

function [fy, dfy] = estimateZfactor(t, y, a, p_pr)   
    fy = -a .* p_pr + (y + y.^2 + y.^3 - y.^4)./(1 - y).^3 - (14.76 .* t - 9.76 .* t.^2 + 4.58 .* t.^3) .* y.^2 + (90.7 .* t - 242.2 .* t.^2 + 42.4 .* t.^3) .*y.^(2.18 + 2.82 .* t);

    dfy = (1 + 4 .* y + 4 .* y.^2   - 4 .* y.^3   + y.^4)./(1 - y).^4 - (29.52 .* t - 19.52 .*   t.^2   + 9.16 .*t.^3)  .* y + (2.18 + 2.82 .* t) .* (90.7 .* t - 242.2 .*   t.^2  + 42.4 .*   t.^3 ) .*   y.^(1.18 + 2.82 .* t);
end




%% Calculates the gas density at pipe conditions
function denGas = gasDensity(sg, temp, pres, zfac)
    Mg = sg.*28.97;                                                % molecular weight of gas
    den_gas_metric  = (pres .* Mg)./(83.143.*zfac.*(temp + 273.15));  % in kg/m.^3
    denGas = 0.0624279606 .* den_gas_metric;                       % in lb/ft.^3 
end

%% Calculates the gas-liquid surface tension
function db = surfaceTension(gas_density, liquid_density)
    db = liquid_density - gas_density;
    db = 15.0 + 0.91.*db;
end


%% gasViscosity: calculates the gas viscosity
function viscGas = gasViscosity(str, t, p, z)
    % E - pipeline   
    % t - temperatures
    % p - pressure
    % z - zfactor
    
    t_r = (t + 273.15.*ones(length(t),1)).*1.8;
    Mg = gasSpecificGravity(str).*28.97;
    
    den_gas = (p.*Mg)./(z.*83.143.*(t + 273.15))./1000;
    
    A1 = ((9.379 + 0.01607.*Mg).*t_r.^1.5)./(209.2  + 19.26.*Mg + t_r);
    A2 = (3.448 + 986.4./t_r + 0.01009.*Mg);
    A3 = 2.447 - 0.2224.*A2;
    
    viscGas = 1e-4.*A1.*exp(A2.*(A2.*den_gas.^A3));
end


%% liquidViscosity: calculates the liquid viscosity
function liqVisc = liquidViscosity(qo, qw, str)
    oil_rate = qo;
    water_rate = qw;
    
    if ((oil_rate + water_rate) < 1e-6)
        oil_rate = 1.0;
    end   
    
    liqVisc = ( oil_rate.*vertcat(str.oil_visc) +water_rate.*vertcat(str.water_visc) )./(oil_rate + water_rate);
end


