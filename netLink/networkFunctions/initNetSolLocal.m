function netSol = initNetSolLocal(wellSol, netSolInit)
    netSolGiven =  (nargin == 2);
    if netSolGiven
        netSol = netSolInit;
    else
        netSol = defaultNetSol(wellSol);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function ns = defaultNetSol                                            %
%%                Creates a default graph for representing the  network.  %
%%                Upward flows have negative sign by convention.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ns = defaultNetSol(ws)    
    ns = struct(...
        'V',   [],...   % set of all vertices
        'Vsrc',   [],...  % set of source vertices in the network
        'Vsnk',   [],... % set of sink vertices in the network
        'Vw',   [],... % set of border nodes connecting with the reservoir
        'VwProd',   [],... % set of production wells
        'VwInj',   [],... % set of injection wells        
        'E',   [], ....  % set of all edges              
        'Eeqp', [], ...  % set of special edges representing equipments in the network                     
        'M',   [], ...   % index matrix mapping source nodes to edges
        'qo',  [], ...   % oil flows in the edges
        'qw',  [], ...   % water flows in the edges
        'qg',  [], ...   % gas flows in the edges
        'pV',  [], ...      % pressures in the nodes
        'boundaryCond', 5*barsa); % network boundary condition (constant pressure at inlet of the separator)
   
    for i = 1:length(ws)        
        if  strcmp(ws(i).name,'')
             error(id('Well:Empty'), ...
                   'Empty well argument is not supported');
        else
            isProducer =   (ws(i).sign == -1);
            isInjector =   (ws(i).sign == 1);            
            if isProducer  % PRODUCER                
                v1 = newVertex(length(ns.V)+1, -1, ws);   
                ns = addVertex(ns,v1, 'isProducer', isProducer);
            elseif isInjector % INJECTOR
                v1 = newVertex(length(ns.V)+1, 1, ws);             
                ns = addVertex(ns,v1, 'isInjector', isInjector);
            else
                error(id('Network:ReservoirConnection'), 'Well interfacing with the reservoir is not a producer nor a injector');
            end
        end
    end   
    
end
