function [netSol] = setWellSolValues(netSol, wellSol, forwardState, p, pScale, varargin)
%SETWELLSOLVALUES set wellSol values in the network
%TODO: handle gas phase flow.

    opt     = struct('ComputePartials',false);                     
    opt     = merge_options(opt, varargin{:});
    
    qWs  = vertcat(wellSol.qWs);
    qOs  = vertcat(wellSol.qOs);        
    pBHP = vertcat(wellSol.bhp);     
    
    pressure  = forwardState.pressure;
    sW = forwardState.s(:,1);
        
 
    if opt.ComputePartials        
%         [~, ~, qWs, qOs, pBHP, p] = initVariablesADI(pressure, sW, qWs, qOs, pBHP, p);
          [~, ~, qWs, qOs, pBHP] = initVariablesADI(pressure, sW, qWs, qOs, pBHP);
    end

    for i=1:length(wellSol)
        well =  getVertex(netSol, netSol.Vw(i));
        well.pressure =  pBHP(i);
        
        well.qoV = qOs(i);
%             well.qgV = wellSol(i).qGs;
        well.qwV = qWs(i);    
        netSol = updateVertex(netSol, well);               
    end
    
    for j=1:numel(netSol.Vc) % controllable vertices
        vertControl = getVertex(netSol, netSol.Vc(j));
        vertControl.pressure = p*pScale;   %% TODO: generalize this using the field 'control' the vertex mock object        
        
        netSol = updateVertex(netSol, vertControl);
        
    end    
    
    %%TODO: update set of controllable edges in createESPNetwork
    for k=1:numel(netSol.Ec) % controllable edges        
       edgeControl = getEdge(netSol, netSol.Ec(k));
       edgeControl.control = p(k);  
       
       netSol = updateEdge(netSol, edgeControl);
    end
   
end

