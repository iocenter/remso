function [netSol] = setWellSolValues(netSol, wellSol, forwardState, varargin)
%SETWELLSOLVALUES set wellSol values in the network
%TODO: handle gas phase flow.

    opt     = struct('ComputePartials',false);                     
    opt     = merge_options(opt, varargin{:});
    
    qWs  = vertcat(wellSol.qWs);
    qOs  = vertcat(wellSol.qOs);        
    pBHP = vertcat(wellSol.bhp);     
    
    p  = forwardState.pressure;
    sW = forwardState.s(:,1);
        
 
    if opt.ComputePartials        
        [~, ~, qWs, qOs, pBHP] = initVariablesADI(p, sW, qWs, qOs, pBHP);
    end

    for i=1:length(wellSol)
        well =  getVertex(netSol, netSol.Vw(i));
        well.pressure =  pBHP(i)*1e-05;        
        
        well.qoV = qOs(i)*day;
%             well.qgV = wellSol(i).qGs;
        well.qwV = qWs(i)*day;    
        netSol = updateVertex(netSol, well);               
    end
end

