function dp =  simpleDp(E, qoE, qwE, qgE, pV)
  
    
    str = vertcat(E.stream);       
    
    %% Producer flows are negative, and injection flows are positive.
    dp = pV*0;
    
    total_rate = qgE + qoE + qwE;   
    flag_rate = abs(double(total_rate)) >= 1e-06*meter^3/day;
    
     if any(~flag_rate)
        if any(flag_rate)
            dp(flag_rate) =  simpleDp(E(flag_rate), qoE(flag_rate), qwE(flag_rate), qgE(flag_rate), pV(flag_rate));
        end
        return 
     end       
 
    
    
    pipes = vertcat(E.pipeline);
    diameters = vertcat(pipes.diam);
    angles = vertcat(pipes.ang);
        

    g = norm(gravity);                         % gravitational constant (ft/s^2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Flow Velocities %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vsl         = superficialLiquidVelocity(qoE, qwE,diameters);           % superficial liquid velocity returns in ft/s    
    vm = vsl;  % in ft/s
   
   %% calculation of the correction factor   
   %TODO: verify next functions
   den_l = liquidDensity(qoE, qwE, str);                          %% liquid density in lb/ft^3           
   
   %% calculating pressure drop due to elevation change
   den_s = den_l;   %% two phase density                          %%  mixture density in lb/ft^3       
   
   dp_el = den_s.* sin(angles).*g;      %% pressure drop due to elevation change  in SI   
%    dp_l = dp_el.*E.pipeline.len;                        %% dp in meters
     
   
   %% calculating friction factor
   den_ns = den_l;      %% no-slip density in lb/ft^3
   vis_ns = liquidViscosity(qoE, qwE, str);      %% no-slip viscosity in SI
   
      
   re_ns = (den_ns.* vm.* diameters)./vis_ns;     %% no-slip reynolds number in SI
   
   roughness = 2.8*10^-5; % in meters   
   friction_factor= (1./ (-1.8 .* log((roughness ./ (diameters) / 3.7).^ 1.11 + 6.9 ./re_ns)./log(10) )).^2;

%  friction_factor =   0.0056 + 0.5./(re_ns).^(0.32);   
   
   dp_f = friction_factor.*den_ns.*vm.^2./(2.*diameters);  % in SI
%    dp_f = dp_f*0;
   
   dp = ((dp_el + dp_f).*(E.pipeline.len)); %% returns in barsa   
   
end

%% superficialLiquidVelocity: calculates the superficial liquid velocity (vsl)
function liqVeloc =  superficialLiquidVelocity(qo, qw, diam)
    
    liquid_rate = (qo + qw);              % liquid rate in sm3/s    
    
    area =  pi.*((diam./2)).^2;                             % pipe radius in m^2
    
    liqVeloc = liquid_rate./area;                 % in m/s
end

%% liqDens: calculates the liquid density (oil and water)
function liqDens = liquidDensity(qoE, qwE, str)
    oil_rate = qoE;
    water_rate = qwE;
   
    if any((oil_rate + water_rate) < 1.e-10)
        warning('Liquid rate approaching to zero. Impossible to calculate liquid density.');
    end    
    oil_mdensity = oil_rate.*vertcat(str.oil_dens);
    water_mdensity =  water_rate.*vertcat(str.water_dens);
    
    
    liqDens = (oil_mdensity + water_mdensity)./(oil_rate + water_rate);   % liquid density in kg/m^3     
end
    
%% liquidViscosity: calculates the liquid viscosity
function liqVisc = liquidViscosity(qo, qw, str)
    oil_rate = qo;
    water_rate = qw;
    
    if ((oil_rate + water_rate) < 1e-8*meter^3/day)
        oil_rate = 1.0;
    end   
    
    liqVisc = ( oil_rate.*vertcat(str.oil_visc) +water_rate.*vertcat(str.water_visc) )./(oil_rate + water_rate);
end