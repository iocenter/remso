function netSol = prodNetwork(wellSol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                        %   
% Creates a simple production network given the installed production and %
% injection wells present in wellSol mock object                         %                                                                       %       
%                                                                        %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    netSol = initNetSolLocal(wellSol);       
    netSol = createSimpleNetwork(netSol);    
end

function netSol = createSimpleNetwork(ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function A = createNetwork: creates the production network by modifying
%                             an adjancency matrix A for the graph 
%                             representing the production network for given
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


