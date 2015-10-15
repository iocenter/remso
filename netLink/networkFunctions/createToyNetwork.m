function netSol = createToyNetwork(ns)
% this functions creates a network containing a very simple production
% infrastructure for each production well.
   %pbar = 14.7*psia;     
    for i=1:length(ns.V)
        isProducer = (ns.V(i).sign == -1);
        isInjector = (ns.V(i).sign == 1);        
        sign = isProducer*-1 + isInjector;
   
        if isProducer
             % well tubing
            pwhV = newVertex(length(ns.V)+1, sign,sign);                       
            ns = addVertex(ns, pwhV);            
            
            prodTubing = newEdge(length(ns.E)+1, ns.V(i), pwhV, sign);
            prodTubing.units = 0; % METRIC=0, FIELD = 1,
            prodTubing.pipeline = newPipeline('diam', 2.5*inch, 'len', 1737*ft , 'ang', degtorad(90), 'temp',  convtemp(175,'F','C'));
            prodTubing.stream = newStream();
            
            prodTubing.stream.gas_visc = 0.0131;
            prodTubing.stream.sg_gas = 0.65;
            prodTubing.stream.gas_dens = 2.6*pound/ft^3;

            prodTubing.stream.oil_visc = 2;
            prodTubing.stream.sg_oil = 0.35;
            prodTubing.stream.oil_dens =  49.9*pound/ft^3;             
            ns = addEdge(ns,prodTubing); 
            
           %% subsea choke as a special equipment            
           
            % horizontal production flowline 
            pdsV = newVertex(length(ns.V)+1, sign,sign);
            ns = addVertex(ns, pdsV);            
            choke = newEdge(length(ns.E)+1, pwhV, pdsV, sign);
%             choke.pipeline = newPipeline('diam', 2.5*inch, 'len', 2737*ft , 'ang', degtorad(90), 'temp',  convtemp(175,'F','C'));
%             choke.stream = prodTubing.stream;
            ns = addEdge(ns, choke, 'isEquipment', true);
                        
            % riser to reach topside facilities
            compV = newVertex(length(ns.V)+1, sign, sign);
            compV.pressure = 30; % 30 bar
            ns = addVertex(ns, compV, 'isSink', true);
        
            prodRiser = newEdge(length(ns.E)+1, pdsV, compV, sign);       
            prodRiser.pipeline = newPipeline('diam', 2.5*inch, 'len', 5737*ft , 'ang', degtorad(90), 'temp',  convtemp(175,'F','C'));
            prodRiser.stream = prodTubing.stream;
            ns = addEdge(ns, prodRiser);
        elseif isInjector  % production well infrastructure      
            %TODO: create injection well infrastructure. Validate pressure
            %drop in the injection pipelines.            
        end         
    end    
    netSol = ns;
end