function [ netSol ] = createSatelliteWellsNetwork(ns)
%createSatelliteWellsNetwork this function creates a really simple network
%for satellite production wells. The script couples two pipes with one
%equipment connecting them to each well.
    numWells = numel(ns.V);    

    sepV = newVertex(length(ns.V)+1, -1);            
    ns = addVertex(ns, sepV, 'isSink', true);

    for i=1:numWells
        isProducer = (ns.V(i).sign == -1);
        isInjector = (ns.V(i).sign == 1);        
        sign = isProducer*-1 + isInjector;
        
        if isProducer
             % well tubing
            pwhV = newVertex(length(ns.V)+1, sign);                       
            ns = addVertex(ns, pwhV);            
            
            prodTubing = newEdge(length(ns.E)+1, ns.V(i), pwhV, sign);
            prodTubing.units = 0; % METRIC=0, FIELD = 1,
            prodTubing.pipeline = wellTubingSettings();
            prodTubing.stream = defaultStream();             
            ns = addEdge(ns,prodTubing); 
            
           %% subsea choke as a special equipment       
            % horizontal production flowline 
            pdsV = newVertex(length(ns.V)+1, sign);
            ns = addVertex(ns, pdsV);            
            choke = newEdge(length(ns.E)+1, pwhV, pdsV, sign);
            ns = addEdge(ns, choke, 'isEquipment', true);
                        
            % riser to reach topside facilities
            prodRiser = newEdge(length(ns.E)+1, pdsV, sepV, sign);       
            prodRiser.pipeline = flowlineRiserSettings();
            prodRiser.stream = defaultStream();
            ns = addEdge(ns, prodRiser);
        end  
    end    
    ns.boundaryCond = 5*barsa; % network boundary condition
    netSol = ns;
end


function [pipe] = wellTubingSettings() %pipeW.dat
     pipe = newPipeline('diam', 76*milli*meter, ... in %m
                      'len', 250, ... % in m
                      'ang', degtorad(90), ...  % in rad
                      'temp', convtemp(60,'C','K'));   % in K  

end

function [pipe] = flowlineRiserSettings() %pipeR.dat
    pipe = newPipeline('diam', 0.24, ... in %m
                      'len', 1000 , ... % in m
                      'ang', degtorad(55), ...  % in rad
                      'temp',  convtemp(60,'C','K'));   % in K  
end

function [str] =  defaultStream() % default stream used in the example
    str = newStream('sg_gas', 0.65, ...  % air = 1
                    'oil_dens', 897,...  % kg/m^3
                    'water_dens', 1025.2, ... % kg/m^3  
                    'oil_visc', 0.00131, ... % Pa s
                    'water_visc', 0.00100); % Pa s
end


function str = gasStream()
    str = newStream();
    str.gas_visc = 0.018*(centi*poise); % viscosity in kilogram/(meter*second)
    str.sg_gas = 0.65;
    str.gas_dens = 2.84*pound/ft^3; % in kg/m^3

    str.oil_visc = 18*(centi*poise); % viscosity in kilogram/(meter*second)
    str.sg_oil = 0.35;    
    str.oil_dens =  56.6*pound/ft^3; % in kg/m^3
end


