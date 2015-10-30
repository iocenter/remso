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
        [~, ~, qWs, qOs, pBHP, p] = initVariablesADI(pressure, sW, qWs, qOs, pBHP, p);
    end

    for i=1:length(wellSol)
        well =  getVertex(netSol, netSol.Vw(i));
        well.pressure =  pBHP(i);
        
        well.qoV = qOs(i);
%             well.qgV = wellSol(i).qGs;
        well.qwV = qWs(i);    
        netSol = updateVertex(netSol, well);               
    end
    
    for j=1:length(netSol.Vc)
        vertControl = getVertex(netSol, netSol.Vc(j));
%         vertControl.pressure = p/barsa;
        vertControl.pressure = p*pScale;
        
        netSol = updateVertex(netSol, vertControl);
        
    end    
   
end

