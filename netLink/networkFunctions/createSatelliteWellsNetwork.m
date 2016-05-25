function [ netSol ] = createSatelliteWellsNetwork(ns)
%createSatelliteWellsNetwork this function creates a really simple network
%for satellite production wells. The script couples two pipes with one
%equipment connecting them to each well.
    for i=1:length(ns.V)
        isProducer = (ns.V(i).sign == -1);
        isInjector = (ns.V(i).sign == 1);        
        sign = isProducer*-1 + isInjector;
   
        if isProducer
             % well tubing
            pwhV = newVertex(length(ns.V)+1, sign);                       
            ns = addVertex(ns, pwhV);            
            
            prodTubing = newEdge(length(ns.E)+1, ns.V(i), pwhV, sign);
            prodTubing.units = 0; % METRIC=0, FIELD = 1,
            prodTubing.pipeline = newPipeline('diam', 0.249*ft, 'len', 1000*meter , 'ang', degtorad(90), 'temp',  convtemp(175,'F','K'));
            prodTubing.stream = gasStream();
            
            prodTubing.stream.gas_visc = 0.0131;
            prodTubing.stream.sg_gas = 0.65;
            prodTubing.stream.gas_dens = 2.6*pound/ft^3;

            prodTubing.stream.oil_visc = 2;
            prodTubing.stream.sg_oil = 0.35;
            prodTubing.stream.oil_dens =  49.9*pound/ft^3;             
            ns = addEdge(ns,prodTubing); 
            
           %% subsea choke as a special equipment            
           
            % horizontal production flowline 
            pdsV = newVertex(length(ns.V)+1, sign);
            ns = addVertex(ns, pdsV);            
            choke = newEdge(length(ns.E)+1, pwhV, pdsV, sign);
            ns = addEdge(ns, choke, 'isEquipment', true);
                        
            % riser to reach topside facilities
            compV = newVertex(length(ns.V)+1, sign);            
            ns = addVertex(ns, compV, 'isSink', true);
        
            prodRiser = newEdge(length(ns.E)+1, pdsV, compV, sign);       
            prodRiser.pipeline = newPipeline('diam', 2.5*inch, 'len', 2000*meter , 'ang', degtorad(120), 'temp',  convtemp(225,'F','K'));
            prodRiser.stream = prodTubing.stream;
            ns = addEdge(ns, prodRiser);
        elseif isInjector  % production well infrastructure      
            %TODO: create injection well infrastructure. Validate pressure
            %drop in the injection pipelines.            
        end         
    end    
    ns.boundaryCond = 5*barsa; % network boundary condition
    netSol = ns;
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


