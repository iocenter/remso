function netSol = runNetwork(ns)
%% runs a full simulation for the whole production network
    
    idsV = ns.Vsrc; % current set of nodes
    Vc = getVertex(ns, idsV);
    
    % Flows from source vertices
    for i=1:length(Vc)
        % flowing from source nodes
        Vin = Vc(i);   
       
        % propagate flows and pressure up to the choke
        [ns, Vin] = propagateFlowPressures(ns, Vin, 'propagPressures', true, 'uptoChoke', true);       

        
        % continue propagating only flows (not pressures) for pipelines  after the chokes
       [ns, Vin] = propagateFlowPressures(ns, Vin, 'propagPressures', false, 'uptoChoke', false);
             
        % from sink node up to downstream the choke back-calculate
       [ns]  = backcalculatePressures(ns,Vin);
    end    
    netSol = ns;
    
    % updating adjacency matrix
%     for i = 1:length(A)  % only the diagonal contains vertices
%         if A(i,i) == -2 %% producer
%             A = runProducer(A,i);                
%         elseif A(i,i) == 2 %% injector
%             A = runInjector(A,i);
%         end        
%     end   
end

function [ns] = backcalculatePressures(ns, Vout, varargin)
% Back-calculate pressures from sink nodes up to downstream the choke.
    Ein = getEdge(ns, Vout.Ein);
    condStop = Ein.equipment;
    while  ~all(condStop)
        % calculating pressure drops in the inlet pipeline                
        [qo, qw, qg, p] = graph2FlowsAndPressures(Vout, Ein);
        dp = dpBeggsBrill(Ein, qo, qw, qg, p);
        Vin = getVertex(ns, Ein.vin);        
        Vin.pressure = Vout.pressure+dp;  % TODO: vector of pressures (equality constraints imposed in the optmizer)
        ns = updateVertex(ns,Vin);
        
        Vout = Vin;
        Ein = getEdge(ns, Vout.Ein);
        if  ~ismember(Vout.id, ns.Vsrc)
            condStop = Ein.equipment;
        else % reached source vertex
            condStop = true;  
        end        
    end
        

end

function [ns, Vin] = propagateFlowPressures(ns, Vin, varargin)
% Propagate flows and pressure up to the choke. Pressure is optional.
    opt     = struct('propagPressures',false, 'uptoChoke', false); % default option    
    opt     = merge_options(opt, varargin{:});


    Eout =  getEdge(ns, Vin.Eout);        
    if opt.uptoChoke
        condStop = Eout.equipment;
    else
        condStop = ismember(Eout.vout, ns.Vsnk);
    end
    
    while  ~all(condStop)
        % updating flows in the edges
        Eout.qoE = Eout.qoE + Vin.qoV;
        Eout.qgE = Eout.qgE + Vin.qgV;
        Eout.qwE = Eout.qwE + Vin.qwV;
        ns = updateEdge(ns,Eout);
        
        % updating flows in the vertices
        Vout = getVertex(ns, vertcat(Eout.vout)); % dest vertices
        Vout.qoV = Vout.qoV + Eout.qoE;
        Vout.qgV = Vout.qgV + Eout.qgE;
        Vout.qwV = Vout.qwV + Eout.qwE;
        ns = updateVertex(ns, Vout);

        % calculating pressure drops in the pipeline        
        if opt.propagPressures           
            [qo, qw, qg, p] = graph2FlowsAndPressures(Vout, Eout);
            dp = dpBeggsBrill(Eout, qo, qw, qg, p);   % TODO: implement BeggsAndBrill for a given inlet pressure.
            Vout.pressure =  Vin.pressure-dp; % TODO: vector of pressures (equality constraints imposed in the optmizer)
            ns = updateVertex(ns,Vout);
        end

        Vin = Vout;
        Eout = getEdge(ns, vertcat(Vin.Eout));
        
        if isempty(Eout)
            condStop = true;
        elseif opt.uptoChoke
            condStop = Eout.equipment;
        end
        
    end    
end

function A = runProducer(A,i)
%% This function performs a flow simulation for a producer
    ik=i;
    vert = A(ik,ik);
    while (vert.type ~= -1) %% production ending vertex
        
        %% updating flows and pressure loss in the pipeline
        colEdge = getOutEdge(A, ik);
        dp = dpin_uphill(A(ik, colEdge), vert);
        A(ik,colEdge).qoE = vert.qoV;
        A(ik,colEdge).qwE = vert.qwV;
        A(ik,colEdge).qgE = vert.qgV;
        A(ik,colEdge).dpE = dp;
    
        %% updating flows and pressure of the ending vertex
        A(colEdge,colEdge).qoV =  A(ik,colEdge).qoE;
        A(colEdge,colEdge).qwV =  A(ik,colEdge).qwE;
        A(colEdge,colEdge).qgV =  A(ik,colEdge).qgE;
        A(colEdge,colEdge).pressure = vert.pressure - dp;
        
        ik = colEdge;
        vert = A(ik,ik);
    end 
end

function A = runInjector(A,j)
%% This function perfoms a simulation for an injection well
    jk = j;
    vert = A(jk,jk);
    while (vert.type ~= 1) %% injection ending vertex
        rowEdge = getInEdge(A, jk);
        dp = dpout_downhill(A(rowEdge, jk), vert);
                
        A(rowEdge,jk).qoE = vert.qoV;
        A(rowEdge,jk).qwE = vert.qwV;
        A(rowEdge,jk).qgE = vert.qgV;
        A(rowEdgeik,jk).dpE = dp;        
        
        
        A(rowEdgeik,rowEdgeik).qoV = vert.qoV;
        A(rowEdgeik,rowEdgeik).qwV= vert.qwV;
        A(rowEdgeik,rowEdgeik).qgV = vert.qgV;
        A(rowEdgeik,rowEdgeik).pressure =  vert.pressure - dp;
        
        jk = rowEdge;
        vert = A(jk,jk);
    end
end