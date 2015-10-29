function dp_psi_tot =  simpleDp(E, qoE, qwE, qgE, pV)
  
    units = struct('METRIC',0, 'FIELD', 1);    
    flow_regime = struct('SEGREGATED', 0, 'TRANSITION', 1, 'INTERMITTENT', 2, 'DISTRIBUTED', 3, 'UNDEFINED', 4);    
    
    str = vertcat(E.stream);       
    
    %% Producer flows are negative, and injection flows are positive.
    dp_psi_tot = pV*0;
    
    total_rate = qgE + qoE + qwE;   
    flag_rate = abs(double(total_rate)) >= 1e-06*meter^3/day;
    
     if any(~flag_rate)
        if any(flag_rate)
            dp_psi_tot(flag_rate) =  simpleDp(E(flag_rate), qoE(flag_rate), qwE(flag_rate), qgE(flag_rate), pV(flag_rate));
        end
        return 
     end       
 
    p_psi       = pV .* (1/psia)*barsa;                                % pressure in psi        
    %p_psi = pV;
    
    
    pipes = vertcat(E.pipeline);
    diameters = vertcat(pipes.diam);
    temperatures = vertcat(pipes.temp);
    angles = vertcat(pipes.ang);
    
    diam_in     = diameters.*(1/inch);                       % pipe diameter in inches  
    

    g = norm(gravity)/(ft/second^2);                         % gravitational constant (ft/s^2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Flow Velocities %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vsl         = superficialLiquidVelocity(qoE, qwE,diam_in);           % superficial liquid velocity returns in ft/s    
    vm = vsl;  % in ft/s
   
   %% calculation of the correction factor   
   %TODO: verify next functions
   den_l = liquidDensity(qoE, qwE, str);                          %% liquid density in lb/ft^3           
   
   %% calculating pressure drop due to elevation change
   den_s = den_l;   %% two phase density                          %%  mixture density in lb/ft^3       
   
   dp_el = den_s.* sin(angles).*g*(ft/second^2);      %% pressure drop due to elevation change  in SI   
%    dp_l = dp_el.*E.pipeline.len;                        %% dp in meters
     
   
   %% calculating friction factor
   den_ns = den_l;      %% no-slip density in lb/ft^3
   vis_ns = liquidViscosity(qoE, qwE, str);      %% no-slip viscosity in SI
   
      
   re_ns = (den_ns.* vm.*(ft/second).* diameters)./vis_ns;     %% no-slip reynolds number in SI
   
   roughness = 2.8*10^-5; % in meters   
   friction_factor= (1./ (-1.8 .* log((roughness ./ (diameters) / 3.7).^ 1.11 + 6.9 ./re_ns)./log(10) )).^2;

%  friction_factor =   0.0056 + 0.5./(re_ns).^(0.32);   
   
   dp_f = friction_factor.*den_ns.*(vm.*(ft/second)).^2./(2.*diameters);  % in SI
%    dp_f = dp_f*0;
   
   dp_psi_tot = ((dp_el + dp_f).*(E.pipeline.len))./barsa; %% returns in barsa   
   
end

%% superficialLiquidVelocity: calculates the superficial liquid velocity (vsl)
function liqRateFt =  superficialLiquidVelocity(qo, qw, diam)
    
    liquid_rate = (qo + qw).*(meter^3/day)./(stb/day);              % liquid rate in bbl/d
    liquid_rate_ft = liquid_rate.*(stb/day)./(ft^3/second);         % liquid rate in ft^3/s    
    
    area =  pi.*((diam./2.*inch/ft)).^2;                             % pipe radius in ft
    
    liqRateFt = liquid_rate_ft./area;                 % in ft/s
end

%% liqDens: calculates the liquid density (oil and water)
function liqDens = liquidDensity(qoE, qwE, str)
    oil_rate = qoE;
    water_rate = qwE;
   
    if any((oil_rate + water_rate) < 1.e-6*meter^3/day)
        warning('Liquid rate approaching to zero. Impossible to calculate liquid density.');
    end
    
    oil_rate_ftd = oil_rate;
    water_rate_ftd = water_rate;
        
    oil_mdensity = oil_rate_ftd.*vertcat(str.oil_dens);
    water_mdensity =  water_rate_ftd.*vertcat(str.water_dens);
    
    
    liqDens = (oil_mdensity + water_mdensity)./(oil_rate_ftd + water_rate_ftd); %returns in SI
%     liqDens =  0.0624279606.*den_metric;                                                        
end
    
%% liquidViscosity: calculates the liquid viscosity
function liqVisc = liquidViscosity(qo, qw, str)
    oil_rate_si = qo*meter^3/day;
    water_rate_si = qw*meter^3/day;
    
    if any((oil_rate_si + water_rate_si) < 1.e-6*meter^3/day)
        warning('Liquid rate approaching to zero. Impossible to calculate liquid density.');
    end
    
    liqVisc = ( oil_rate_si.*vertcat(str.oil_visc) +water_rate_si.*vertcat(str.water_visc) )./(oil_rate_si + water_rate_si); % in kg/(m*s)
    
end