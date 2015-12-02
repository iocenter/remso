function netSol = runNetwork(ns, wellSol, forwardState,p, pScale,  varargin)
%% runs a full simulation for the whole production network
    
    opt     = struct('ComputePartials',false);                     
    opt     = merge_options(opt, varargin{:});

    ns = setWellSolValues(ns, wellSol, forwardState, p, pScale, 'ComputePartials',opt.ComputePartials);

    idsV = ns.Vsrc; % current set of nodes
    Vsrc = getVertex(ns, idsV);
    
    % Flows from source vertices
    for i=1:length(Vsrc)
        % flowing from source nodes
        Vin = Vsrc(i);   
       
        % propagate flows and pressure up to the choke
        [ns, Vin] = propagateFlowPressures(ns, Vin, 'propagPressures', true, 'uptoChokeOrPump', true);       

        
        % continue propagating only flows (not pressures) for pipelines  after the chokes
       [ns, Vin] = propagateFlowPressures(ns, Vin, 'propagPressures', false, 'uptoChokeOrPump', false);      
    end 
    
    surfaceSinks = setdiff(vertcat(ns.Vsnk), vertcat(ns.VwInj));
    for j=1:numel(surfaceSinks)      
        Vout = getVertex(ns, surfaceSinks(j));
        
        % from sink node up to downstream the choke back-calculate
        [ns]  = backcalculatePressures(ns,Vout);
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

function pin = dpPressurePipes(Ein, Vout)
    % calculating pressure drops in the inlet pipeline                
    [qo, qw, qg, p] = graph2FlowsAndPressures(Vout, Ein);    
    voutP = vertcat(Vout.pressure);
    %%TODO: be sure that dp and pin are ADIs if necessary.
    dp = qo*0;
    pin = qo*0;
    
    for i=1:numel(Ein)
        condEsp = (vertcat(Ein.esp));    
        if any(condEsp)
            str = Ein(i).stream;
            
            totalFlow = vertcat(Ein(i).qoE) + vertcat(Ein(i).qwE);
            avgDen = (str.water_dens + str.oil_dens)/2;
            dp(i) = pump_dp_cst(totalFlow, avgDen);
            
            %% TODO: consider the impact of the frequency in the dp of the pump
            pin(i) = voutP(i) + dp(i);
        end
        
        if any(~condEsp)
            dp(i) = dpBeggsBrill(Ein(i), qo(i), qw(i), qg(i), p(i));                           
            pin(i) = voutP(i)+dp(i);
        end
    end
    
end

function pin = dpPressureEquip(Ein, Vout)
    pin = vertcat(Vout.pressure);
%     condEsp = (vertcat(Ein.esp));    
%     if any(condEsp)
%         str = Ein(condEsp).stream;
%         
%         totalFlow = vertcat(Ein(condEsp).qoE) + vertcat(Ein(condEsp).qwE);
%         avgDen = (str.water_dens + str.oil_dens)/2;
%         dp = pump_dp_cst(totalFlow, avgDen);
%         
%         %% TODO: consider the impact of the frequency in the dp of the pump
%         pin(condEsp) = vertcat(Vout(condEsp).pressure) + dp;
%     end
end
    
function newV = updatePressures(v, pres)
    for i=1:numel(v)
        v(i).pressure = pres(i);
    end
    newV = v;
end

function [ns] = backcalculatePressures(ns, Vout, varargin)
% Back-calculate pressures from sink nodes up to downstream the choke.
    Ein = getEdge(ns, Vout.Ein);
    vin = getVertex(ns, Ein.vin);
    
    condStop = vin.flagStop;
    while  ~all(condStop)        
        condEquip = (vertcat(Ein.equipment));
         if numel(Vout) == 1 &&  numel(Ein) > 1 %% it is a manifold
            Vout = repmat(Vout, numel(condEquip), 1);
        end        
        if any(condEquip)
              %% TODO: error here !! updating inlet pressure of the pipe (equip) with the outlet pressure!!  
%             pres = dpPressureEquip(Ein(condEquip), Vout(condEquip));
%             vin(condEquip) = updatePressures(vin(condEquip), pres);

%             if Ein.separator  %% TODO: remove part of the water according to the separator efficiency               
%             end
  
            
        elseif any(~condEquip)                    
            pres = dpPressurePipes(Ein(~condEquip), Vout(~condEquip));
            vin(~condEquip) = updatePressures(vin(~condEquip), pres);
            
        end
        
%         if Ein.equipment            
%             vin.pressure = Vout.pressure;             
%         else
%             % calculating pressure drops in the inlet pipeline                
%             [qo, qw, qg, p] = graph2FlowsAndPressures(Vout, Ein);
%             dp = dpBeggsBrill(Ein, qo, qw, qg, p);               
%             vin.pressure = Vout.pressure+dp;                                    
%         end
        ns = updateVertex(ns,vin);  
        Vout = vin;
        Ein = getEdge(ns, vertcat(Vout.Ein));
        vin = getVertex(ns, vertcat(Ein.vin));
        condStop = vertcat(vin.flagStop);         
    end
        

end

function [ns, Vin] = propagateFlowPressures(ns, Vin, varargin)
% Propagate flows and pressure up to the choke. Pressure is optional.
    opt     = struct('propagPressures',false, 'uptoChokeOrPump', false); % default option    
    opt     = merge_options(opt, varargin{:});


    Eout =  getEdge(ns, Vin.Eout);        
    if opt.uptoChokeOrPump
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
            [qo, qw, qg, p] = graph2FlowsAndPressures(Vin, Eout);
            dp = dpBeggsBrill(Eout, qo, qw, qg, p);
            Vout.pressure =  Vin.pressure-dp;
            Vout.flagStop = true;
            ns = updateVertex(ns,Vout);
        end

        Vin = Vout;
        Eout = getEdge(ns, vertcat(Vin.Eout));
        
        if isempty(Eout)
            condStop = ones(length(Eout),1);
        elseif opt.uptoChokeOrPump
            condStop = Eout.choke || Eout.pump;
        else
            condStop = zeros(length(Eout),1);
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