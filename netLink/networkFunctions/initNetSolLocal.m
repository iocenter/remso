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
        'Vint',   [],... % set of interior vertices in the network
        'Vc',   [],... % set of control vertices in the network
        'E',   [], ....  % set of all edges      
        'Ec',  [], ...   % set of controllable edges
        'Eeqp', [], ...  % set of special edges representing equipments in the network                                 
        'Echk', [], ...  % subset of special edges denoting the chokes of the network
        'Epmp', [], ...  % subset of special edges denoting the pumps of the network.
        'Esep', [], ...  % subset of special edges denoting separator in the network.
        'Esrc', [], ...  % set of edges leaving a source node in Vsrc
        'Esnk', [], ...  % set of edges reaching a sink node in Vsnk
        'A',   []);      % incidency matrix 
        
    for i = 1:length(ws)        
        if  strcmp(ws(i).name,'')
             error(id('Well:Empty'), ...
                   'Empty well argument is not supported');
        else
            isProducer =   (ws(i).sign == -1);
            isInjector =   (ws(i).sign == 1);            
            if isProducer  % PRODUCER                
                v1 = newVertex(length(ns.V)+1, -1, -2, ws);   
                ns = addVertex(ns,v1, 'isProducer', isProducer);
            elseif isInjector % INJECTOR
                v1 = newVertex(length(ns.V)+1, 1, 2, ws);             
                ns = addVertex(ns,v1, 'isInjector', isInjector);
            else
                error(id('Network:ReservoirConnection'), 'Well interfacing with the reservoir is not a producer nor a injector');
            end
        end
    end
    
    ns.A = initAdjancencyMatrix(ns);
end

function A = initAdjancencyMatrix(ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function A = initAdjancencyMatrix(ns): init the adjacency matrix with  %
%%                                        producers and injectors.        %    
%%                                                                        %     
%%          Main diagonal of A: -2 represents a producer,                 %  
%%                               2 represents an injector.                %  
%%                              -1 and 1 represent final vertices that are% 
%%                                not in connection with reservoir.       %
%%                                                                        %  
%%          Rows: represents the nodes sending production.                %
%%               Ex: an edge connecting v1 to v2 is positioned in A(v1,v2)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = zeros(length(ns.V));
    for i=1:length(ns.V)
        A(i,i) = ns.V(i).id;
    end
end

function A = createAdjacencyMatrixFromFile(vk, netinput)
%%TODO: finish this implementation using eclipse input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function A = createAdjacencyMatrix: creates an adjancency matrix A for % 
%%                                     the graph representing the prod    %
%%                                     nettwork algebraic variables vk.   %
%%                                                                        %     
%%          Main diagonal of A: -2 represents a producer,                 %  
%%                               2 represents an injector.                %  
%%                              -1 and 1 represent final vertices that are% 
%%                                not in connection with reservoir.       %
%%                                                                        %  
%%          Rows: represents the nodes sending production.                %
%%               Ex: an edge connecting v1 to v2 is positioned in A(v1,v2)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %% vertices columns: vert    id  sign    type    press
    %% edges columns: edge    vo  vd  sign    diam    len ang temp

    fid = fopen('networkInput.data');
    tline = fgets(fid);
    if ~strcmp(tline,'$')
        return;
    end

    nV = str2double(fgets(fid));
    A = zeros(nv);
    for i = 1:nV
        tline = fgets(fid);
        str = strsplit(tline, ' ');

        numVert = str2double(str(1));
        idVert = str(2);
        signVert = str(3);
        typeVert  = str(4);
        pressVert = str(5);

        A(numVert, numVert) = typeVert;        

        %% TODO: update vertex information

    end
    nE = str2double(fgets(fid));

    for j=1:nE
        tline = fgets(fid);
        str = strsplit(tline, ' ');


        nameEdge = str(1);
        vorigin = str(2);
        vdestiny = str(3);
        sign = str(4);

        diam = str(5);
        lengt = str(6);
        angle = str(7);
        temp  = str(8);

        %% TODO: create edge mock object


        %% TODO: insert it into the adjacency matrix

    end
    tline = fgets(fid);
    if ~strcmp(tline,'$')
        return; %error
    end

    fclose(fid);
end

