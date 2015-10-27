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



function dp =  dpBeggsBrill(E, qoE, qwE, qgE, pV)
%% dp: calculates the pressure drop for inlet flow in psi
%% TODO: extend this function to a vector of edges E = [e1, e2,.. e_m]
     
%     units = struct('METRIC',0, 'FIELD', 1);    
    flow_regime = struct('SEGREGATED', 0, 'TRANSITION', 1, 'INTERMITTENT', 2, 'DISTRIBUTED', 3, 'UNDEFINED', 4);    
    
    str = vertcat(E.stream);       
    
    %% Producer flows are negative, and injection flows are positive.
    dp = pV*0;
    zfactor = pV*0;
    
    total_rate = qgE + qoE + qwE;   
    flag_rate = abs(double(total_rate)) >= 1e-06*meter^3/day;
    flag_gas = abs(double(qgE)) >= 1e-05*meter^3/day;
%     assert(~any(flag_gas));
    
    
     if any(~flag_rate)
        if any(flag_rate)
            dp(flag_rate) =  dpBeggsBrill(E(flag_rate), qoE(flag_rate), qwE(flag_rate), qgE(flag_rate), pV(flag_rate));
        end
        return 
     end       
 
    
    pipes = vertcat(E.pipeline);
    diameters = vertcat(pipes.diam);
    temperatures = vertcat(pipes.temp);
    angles = vertcat(pipes.ang);        
                               
    grav = norm(gravity);                                       % gravitational constant (m/s^2)
    if any(flag_gas)            
        strGas = vertcat(str.sg_gas);
        
        zfactor(flag_gas) = gasZFactor(strGas(flag_gas), temperatures(flag_gas), pV(flag_gas));     % gas z-factor
    end   
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Flow Velocities %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vsl         = superficialLiquidVelocity(qoE, qwE,diameters);           % superficial liquid velocity returns in ft/s
    if any(flag_gas)
        assert(numel(flag_gas)==1);    
        vsg         = superficialGasVelocity(qgE, pV, zfactor, diameters, temperatures); % superficial two phase velocity in ft/s
        vm          = vsl + vsg;                                     % superficial two phase velocity in ft/s
    else
        vm = vsl;
        vsg=0;
    end
    
    froude_num = vm.^2./diameters./grav;                                 % froude number    
    liquid_content = vsl./vm;
    
    
    %% if the liquid content is 0, changing it to something small
%     liquid_content = max(liquid_content,  1e-8.*ones(length(E),1));
    if any(liquid_content < 1e-8.*ones(length(E),1))
        warning('Liquid flow is approaching to zero.');
    end
        
        
    %% constants used in formula for evaluating the flow regime  
    l1 = 316.*liquid_content.^(0.302);
    l2 = 0.0009252.*liquid_content.^(-2.4684);
    l3 = 0.10.*liquid_content.^(-1.4516);
    l4 = (0.5).*(liquid_content.^(-6.738));
    
    %%%%% BaBr.pdf report %%%
%     l1 = 230;
%     l2 = 0.0124;
%     l3 = 0.456;
%     l4 = 590;
%     liquid_content = 0.35;
%     froude_num = 29.6;
    %%%%% BaBr.pdf report %%%
        
      
    %% conditional vectors which are used to determine flow regimes
    cond_seg = (liquid_content < 0.01 & froude_num < l1) | (liquid_content>=0.01 & froude_num < l2);
    cond_trans = (liquid_content >= 0.01 & froude_num >= l2 & froude_num <= l3);    
    cond_dist = (liquid_content < 0.4 & froude_num >= l1) | (liquid_content >= 0.4 & froude_num > l4);
    cond_interm = (liquid_content >= 0.01 & liquid_content < 0.4 & froude_num > l3 & froude_num <= l1) | ...
                  (liquid_content >= 0.4 & froude_num > l3 & froude_num <= l4) | ...
                  (~cond_seg & ~cond_trans & ~cond_dist);  % assuming intermittent flow when it is not possible to determine.

   
   
   
   %% calculation of the correction factor      
   den_l = liquidDensity(qoE, qwE, str);                          %% liquid density in SI
   den_g = pV*0;
   if  any(flag_gas)       
      den_g(flag_gas) = gasDensity(gasSpecificGravity(str(flag_gas)), temperatures(flag_gas), pV(flag_gas), zfactor(flag_gas));    %% gas density in kg/m^3
%       den_g(flag_gas) = 2.6; % Babr.pdf example
   end
   surface_tens = surfaceTension(den_g, den_l);          %% gas - liquid surface tension in SI
   
   nlv = vsl .* (den_l./(grav.*surface_tens)).^(0.25);         %% liquid velocity number
   
   %% Payne correction factor to holdup  (Page 43 of the Book 'James Brill, Multiphase Flow in Wells, No 17')                               
   payne_cor = 0.924.*ones(length(E),1);
   
   % Beggs % Brill holdup constants
   regime = zeros(length(E),1); 
   cor    = zeros(length(E),1);
   
   a = zeros(length(E),1);
   b = zeros(length(E),1);
   c = zeros(length(E),1);
   
   d = zeros(length(E),1);
   e = zeros(length(E),1);
   f = zeros(length(E),1);
   g = zeros(length(E),1);  
   
   condAng = (angles < zeros(length(angles),1));
   if any(condAng) % Negative Angle
       a(condAng) = 0.98;
       b(condAng) = 0.4846;
       c(condAng) = 0.0868;
       
       d(condAng) = 4.7;
       e(condAng) = -0.3692;
       f(condAng) = 0.1244;
       g(condAng) = -0.5056;
  
       error('TODO: correct ABC constants for downhill flow regimes.')
       payne_cor(condAng) = 0.685;
   end   
   
   if any(~condAng) && any(cond_seg) % SEGREGATED Flow        
       a(cond_seg) = 0.98;
       b(cond_seg) = 0.4846;
       c(cond_seg) = 0.0868;
       
       d(cond_seg) = 0.011;
       e(cond_seg) = -3.768;
       f(cond_seg) = 3.539;
       g(cond_seg) = -1.614;        
       
       regime(cond_seg) =  flow_regime.SEGREGATED;
   end
   
   if any(~condAng) &&any(cond_interm) % INTERMITTENT flow
       
       a(cond_interm) = 0.845;
       b(cond_interm) = 0.5351;
       c(cond_interm) = 0.0173;
       
       d(cond_interm) = 2.96;
       e(cond_interm) = 0.305;
       f(cond_interm) = -0.4473;
       g(cond_interm) = 0.0978;
       
       regime(cond_interm) = flow_regime.INTERMITTENT;              
  
   end
   if any(~condAng) && any(cond_dist) % DISTRIBUTED flow
       
       a(cond_dist) = 1.065;
       b(cond_dist) = 0.5824;
       c(cond_dist) = 0.0609;
       
       d(cond_dist) = 0;
       e(cond_dist) = 0;
       f(cond_dist) = 0;
       g(cond_dist) = 0;
       
       regime(cond_dist) = flow_regime.DISTRIBUTED;
       
   end
   
   if any(~condAng) && (cond_trans) % TRANSITION flow
       regime(cond_trans) = flow_regime.TRANSITION;   
       error('Transition flow is not handled.')
   end      
   
   if any(~condAng) && any(~cond_dist)
       cor(~cond_dist) = (1 - liquid_content(~cond_dist)).*log(d.*(liquid_content(~cond_dist).^e).*(nlv(~cond_dist).^f).*(froude_num(~cond_dist).^g));
   end
   
   %% Horizontal liquid holdup, hl(0), except for TRANSITION flow,
   %% which will be later calculated as an interpolation of SEGREGATED and
   %% INTERMITTENT flows.   
   hz_holdup = (a.*(liquid_content.^b))./(froude_num.^c);

   assert(all(or(cond_interm,or(cond_seg,or(cond_trans,cond_dist)))),'Couldn''t find flow regime');   
   
    %% the horizontal holdup is equal to liquid content if smaller
   hz_holdup = max(hz_holdup, liquid_content);     

   %% checking if correction >= 0
   condCor = (cor < 0);
   if any(condCor)
%        disp '### Warning ###'
%        disp 'From: Beggs & Brill 1973'
%        disp 'The calculated correction factor, C, is negative... '
%        disp 'Resetting to 0.0'
       warning('The calculated correction factor, C, is negative... Resetting to 0.0');
       cor(condCor) = zeros(sum(condCor),1);
   end
   
   %% liquid holdup correction for inclination
   if any(~condAng) && any(~cond_dist)       
       psi(~cond_dist) = ones(length(angles(~cond_dist)),1) + cor(~cond_dist).*(sin(1.8.*angles(~cond_dist)) - 0.333.*(sin(1.8.*angles(~cond_dist)).^3));
   elseif any(~condAng) && any(cond_dist)
        psi(cond_dist) = 1;
   end
       
   holdup = payne_cor.*hz_holdup.*psi; % liquid holdup corrected for inclination
   
   %% if trasition regime, the liquid holdup is a mix of segregated and intermmitent
   if any(cond_trans)
       holdup(cond_trans) = liqholdupTransitionFlow(liquid_content(cond_trans), nlv(cond_trans), froude_num(cond_trans), angles(cond_trans), l2(cond_trans), l3(cond_trans));
   end
   
   %% corrects the holdup if wrong
   condHoldup = holdup > 1;
   if any(condHoldup)        
      holdup = -(max(-holdup, -1.0.*ones(length(holdup),1)));
      warning('Holdup is greater than 1.');
   end      
   
   %% two phase density
   den_s = den_l .* holdup + den_g .* (1 - holdup);  
   
   %% calculating pressure drop due to elevation change
   dp_el = den_s.*norm(grav).*sin(angles);    %% pressure gradient due to elevation change

   %% calculating friction factor
   vis_g = 0*pV;
   if any(flag_gas)
       vis_g(flag_gas) = gasViscosity(str(flag_gas), temperatures(flag_gas), pV(flag_gas), zfactor(flag_gas));
   end
       
   den_ns = den_l .* liquid_content + den_g .* (1 - liquid_content);      %% no-slip density
   vis_ns = liquidViscosity(qoE, qwE, str) .* liquid_content + vis_g .* (1 - liquid_content);       %% no-slip viscosity   
   
   re_ns = (den_ns .* vm .* diameters)./vis_ns;     %% no-slip reynolds number
   
   
   %%TODO: why denominator is log(-x) ? This is producing an irrational
   %%number.
   % reynolds_threshold = 10^(3.8215/4.5223) ~≃ 7
   reynolds_threshold = 7;
   condR = (re_ns > reynolds_threshold);  %% codas' number
   fn = re_ns;
   if any(condR)
       fn(condR) = 1./(2 .* log((re_ns(condR) ./ (4.5223 .* log(re_ns(condR))./log(10) - 3.8215)))./log(10)).^2;    %% no-slip friction factor
   end
   if any(~condR)
       fn(~condR) = 0.0056 + 0.5./(re_ns(~condR)).^(0.32); % simple calculation for fn
   end
   
   cond_fn0 = fn < 0;
   if any(cond_fn0)
       warning('Friction factor is negative. Setting up to zero.')
       fn(cond_fn0) = 0;       
   end
   y = liquid_content ./(holdup.^2);
   
   cond_y = (y > 1.0 & y < 1.2);   
   s_term = y; % initialization, does not require to be y.
   if any(cond_y)
       s_term(cond_y) = log(2.2.*y(cond_y) -1.2); 
   end
   if any(~cond_y)
       s_term(~cond_y) = log(y(~cond_y)) ./ (-0.0523 + 3.182 .* log(y(~cond_y)) - 0.8725 .* (log(y(~cond_y)).^2) + 0.01853 .* (log(y(~cond_y)).^4));
   end

   ftp = fn .* exp(s_term);      %% the friction factor

   %% calculating pressure drop due to friction   
   dp_f = (ftp .* den_ns .* (vm.^2)) ./ (2.*diameters);


   %% calculating acceleration term
   ek = 2.16e-4 .* (den_ns .* vm .* vsg) ./ pV;  %% it has impact only when there is gas flowing in the pipe

   %% calculating total pressure drop
   
   dp_tot = (dp_f + dp_el) ./ (1 - ek);  %% total pressure drop per length of pipe (in Pa./ m)
   %% double dp_tot_bar = dp_tot ./ 14.5038 ./ 0.3048;       // total pressure drop in bar / m


   %% converting length from m
   lengths =  vertcat(pipes.len);  
   

   %% total pressure drop in Pascal
   dp = dp_tot .* lengths;    

end

function [holdup] = liqholdupTransitionFlow(liquid_content, nlv, froude_num, angles, l2, l3)
       frac = (l3-froude_num)./(l3-l2);
       
       %% SEGREGATED flow constants
       a_seg = 0.98;
       b_seg = 0.4846;
       c_seg = 0.0868;
       
       %% INTERMITTENT flow constants
       a_int = 0.845;
       b_int = 0.5351;
       c_int = 0.0173;
       
       %% horizontal holdups
       hz_holdup_seg = frac.*(a_seg.*(liquid_content.^b_seg))./(froude_num.^c_seg);       
       hz_holdup_int = (1-frac).*(a_int.*(liquid_content.^b_int))./(froude_num.^c_int);
       
       %% horizontal holdup is the liquid content if smaller.
       %% -(max(-x,-y)) = min(x,y) for positive integers
       hz_holdup_seg = -(max(-hz_holdup_seg, -liquid_content));       
       hz_holdup_int = -(max(-hz_holdup_int, -liquid_content));
       
       %% correction factors
       cor_seg = zeros(length(liquid_content),1);       
       cor_int = zeros(length(liquid_content),1);
       
       cond_ang = (angles < 0);
       if any(cond_ang)           
           d_ang = 4.7;
           e_ang = -0.3692;
           f_ang =  0.1244;
           g_ang = -0.5056;           
           
           cor_seg(cond_ang) = (1-liquid_content(cond_ang)).*log(d_ang.*(liquid_content(cond_ang).^e_ang).*(nlv(cond_ang).^f_ang).*(froude_num(cond_ang).^g_ang));
           cor_int(cond_ang) = cor_seg(cond_ang);
       elseif any(~cond_ang)           
           d_seg = 0.011;
           e_seg = -3.768;
           f_seg = 3.539;
           g_seg = -1.614;
           
           d_int = 2.96;
           e_int = 0.305;
           f_int = -0.4473;
           g_int = 0.0978;
           
           cor_seg(~cond_ang) = (1 - liquid_content(~cond_ang)) .* log(d_seg.*(liquid_content(~cond_ang).^e_seg).*(nlv(~cond_ang).^f_seg).*(froude_num(~cond_ang).^g_seg));
           cor_int(~cond_ang) = (1 - liquid_content(~cond_ang)) .* log(d_int .* (liquid_content(~cond_ang).^e_int) .* (nlv(~cond_ang).^f_int) .*(froude_num(~cond_ang).^g_int));
       end 
       
       psi_seg = 1 + cor_seg .* (sin(1.8.*angles) - 0.333 .* (sin(1.8.*angles)).^3 );
       psi_int = 1 + cor_int .* (sin(1.8.*angles) - 0.333 .* (sin(1.8.*angles)).^3);

%        holdup = payne_cor .* (frac .* (hz_holdup_seg .* phi_seg) + (1 - frac) .* (hz_holdup_int .* phi_int));
       holdup = (frac .* (hz_holdup_seg .* psi_seg) + (1 - frac) .* (hz_holdup_int .* psi_int));
end

%% gasSpecificGravity: returns the gas specific gravity
function sg = gasSpecificGravity(str)
    sg = vertcat(str.sg_gas);
end

%% superficialGasVelocity: calculates the superficial gas velocity (vsg)
function gv = superficialGasVelocity(qgE, pres, zfac, diam, temp)    
    gas_rate_surface = qgE/(ft^3/second);  % the gas rate in Sft^3/s
    
    % surface conditions
    Tsc = convtemp(60,'F','R');
    Psc = 14.7; % in psia
    
    % pipe conditions
    temp_R = convtemp(temp,'F','R');
                                                                                                                    
    A = pi.*((diam./2).*inch/ft).^2;     % pipe area in ft^2                                                                     
    
    gv = (gas_rate_surface.*zfac.*temp_R.*Psc)./(A.*Tsc.*pres); % in ft/s    
    gV = gv*(ft/second); % in m/s
end


%% superficialLiquidVelocity: calculates the superficial liquid velocity (vsl)
function liqVeloc =  superficialLiquidVelocity(qo, qw, diam)
    
    liquid_rate = (qo + qw);              % liquid rate in sm3/s    
    
    area =  pi.*((diam./2)).^2;                             % pipe radius in m
    
    liqVeloc = liquid_rate./area;                 % in m/s
end

%% liqDens: calculates the liquid density (oil and water)
function liqDens = liquidDensity(qoE, qwE, str)
    oil_rate = qoE;
    water_rate = qwE;
   
    if any((oil_rate + water_rate) < 1.e-8)
        warning('Liquid rate approaching to zero. Impossible to calculate liquid density.');
    end    
    oil_mdensity = oil_rate.*vertcat(str.oil_dens);
    water_mdensity =  water_rate.*vertcat(str.water_dens);
    
    
    liqDens = (oil_mdensity + water_mdensity)./(oil_rate + water_rate);   % liquid density in kg/m^3     
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
    t = convtemp(t,'K','F');                %% K to F
    p = p./barsa;                           %% Pascal to bar
    
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
    air_molecular_weight = 28.97*gram/(kilogram);                                    % in kg/mol
    ideal_gas_cons = 8.3144598;                                                      % J/mol K

    Mg = sg.*air_molecular_weight;                                                   % molecular weight of gas
    
    denGas  = (pres .* Mg)./(ideal_gas_cons.*zfac.*temp);                    % in kg/m.^3
end

%% Calculates the gas-liquid surface tension
function db = surfaceTension(gas_density, liquid_density)
    %% convertign to original unity in resopt to perform calculations
    
    gas_density_metric = gas_density./(pound/(ft^3));               % gas density in lb/ft^3
    liq_density_metric = liquid_density./(pound/(ft^3));            % liq density in lb/ft^3


    db_metric = liq_density_metric - gas_density_metric;
    db_metric = 15.0 + 0.91.*db_metric;                             % surface tension in dyn/cm
    
    db = db_metric.*(dyne)./(centi*meter);                          % surface tension in J / m    
end


%% gasViscosity: calculates the gas viscosity
function viscGas = gasViscosity(str, t, p, z)
    % E - pipeline   
    % t - temperatures
    % p - pressure
    % z - zfactor
    
    t_r = t.*1.8;
    R = 8.3143;    
    molecular_weight_air = 28.97.*gram;
    Mg = gasSpecificGravity(str).*molecular_weight_air;    
    
    den_gas = (p.*Mg)./(z.*R.*(t));
    
    A1 = ((9.379 + 0.01607.*Mg).*t_r.^1.5)./(209.2  + 19.26.*Mg + t_r);
    A2 = (3.448 + 986.4./t_r + 0.01009.*Mg);
    A3 = 2.447 - 0.2224.*A2;
    
    viscGas = 1e-4.*A1.*exp(A2.*(A2.*den_gas.^A3)); %% TODO: error in the exponential formulae !!!
end


%% liquidViscosity: calculates the liquid viscosity 
function liqVisc = liquidViscosity(qo, qw, str)
    oil_rate = qo;
    water_rate = qw;
    
    if ((oil_rate + water_rate) < 1e-8)
        oil_rate = 1.0;
    end   
    
    liqVisc = ( oil_rate.*vertcat(str.oil_visc) +water_rate.*vertcat(str.water_visc) )./(oil_rate + water_rate);
end


