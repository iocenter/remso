function netSol = prodNetwork(wellSol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                        %   
% Creates a simple production network given the installed production and %
% injection wells present in wellSol mock object                         %                                                                       %       
%                                                                        %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    netSol = initNetSolLocal(wellSol);       
    %netSol = createSimpleNetwork(netSol);    
    
    netSol = createToyNetwork(netSol);
end

function netSol = createToyNetwork(ns)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Edge e1 %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            compV.pressure = 800;
            ns = addVertex(ns, compV, 'isSink', true);
        
            prodRiser = newEdge(length(ns.E)+1, pdsV, compV, sign);       
            prodRiser.pipeline = newPipeline('diam', 2.5*inch, 'len', 2737*ft , 'ang', degtorad(90), 'temp',  convtemp(175,'F','C'));
            prodRiser.stream = prodTubing.stream;
            ns = addEdge(ns, prodRiser);
        elseif isInjector  % production well infrastructure      
            %TODO: create injection well infrastructure. Validate pressure
            %drop in the injection pipelines.            
        end         
    end    
    netSol = ns;
end


function netSol = createSimpleNetwork(ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function A = createNetwork: creates the production network by modifying
%                             an adjancency matrix A for the graph 
%                              representing the production network for given
%                             algebraic variables nk
%
%               Main diagonal of A: -2 represents a producer, 2 represents
%                                   an injector.
%                                   -1 and 1 represent final vertices 
%                                   that are not in connection with reservoir.
%
%               Rows: represents the nodes sending production. 
%               Ex: an edge connecting v1 to v2 is positioned in A(v1,v2).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    for i=1:length(ns.V)
        isProducer = (ns.V(i).sign == -1);
        isInjector = (ns.V(i).sign == 1);        
        sign = isProducer*-1 + isInjector;
        
        % basic injection well infrastructure                
        if isInjector 
            % tubing from wellbore up to the  wellhead
            whV = newVertex(length(ns.V)+1, sign,  sign);
            ns = addVertex(ns, whV);
            
            %% TODO: correct tubing angle (loop at pipesim model)
            tubing = newEdge(length(ns.E)+1, whV, ns.V(i), sign);
            tubing.pipeline = newPipeline('diam', 139.7*(10^3) , 'len',  2590.8,'ang',0*(pi/180),  'temp', 4.44*(274.15));
            ns = addEdge(ns,tubing);
            
             % horizontal flowline
            conV = newVertex(length(ns.V), sign, sign);
            ns = addVertex(ns, conV);
             
            flowline = newEdge(length(ns.E)+1, conV, whV, sign );
            flowline.pipeline =  newPipeline('diam',101.6*(10^3),  'len', 365.76, 'ang', 90*(pi/180), 'temp', 7.2*(274.15));
            ns = addEdge(ns, flowline);
            
            % riser after the topside compressor        
            compV = newVertex(length(ns.V)+1, sign, sign);
            ns = addVertex(ns, compV);
            
            riser = newEdge(length(ns.E)+1, compV, whV, sign);
            riser.pipeline = newPipeline('diam', 101.6*(10^3), 'len', 76.2, 'ang', 0*(pi/180), 'temp', 7.2*(274.15));
            ns = addEdge(ns, riser);
            
        elseif isProducer  % production well infrastructure       
            % well tubing
            pwhV = newVertex(length(ns.V)+1, sign,sign);            
            ns = addVertex(ns, pwhV);
            
            prodTubing = newEdge(length(ns.E)+1, ns.V(i), pwhV, sign);
            prodTubing.pipeline = newPipeline('diam', 139.7*(10^3), 'len', 3535 , 'ang', 35.6*(pi/180), 'temp',  15.56*(274.15));
            ns = addEdge(ns,prodTubing);
            
            % horizontal production flowline 
            prodconV = newVertex(length(ns.V)+1, sign,sign);
            ns = addVertex(ns, prodconV);
            
            prodFlowline = newEdge(length(ns.E)+1, pwhV, prodconV, sign);
            prodFlowline.pipeline = newPipeline('diam', 131.7498*(10^3), 'len', 3000 , 'ang', 90*(pi/180), 'temp',  15.56*(274.15));
            ns = addEdge(ns, prodFlowline);
            
            
            % riser to reach topside facilities
            compV = newVertex(length(ns.V)+1, sign, sign);
            ns = addVertex(ns, compV);
        
            prodRiser = newEdge(length(ns.E)+1, ns.V(i), compV, sign);       
            prodRiser.pipeline = newPipeline('diam', 131.7498*(10^3), 'len', 24.384 , 'ang', 0*(pi/180), 'temp',  15.56*(274.15));
            ns = addEdge(ns, prodRiser);
        else
            error('Error: Vertex interfacing reservoir and network is not a well !');        
        end        
       
    end
    
    netSol = ns;
end


