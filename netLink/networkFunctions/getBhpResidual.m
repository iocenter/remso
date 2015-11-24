function [ dp] = getBhpResidual(wellSol, netSol)
%GETBHPRESIDUAL get bottom-hole pressure residual (difference from reservoir and network simulations)
    Epmp = getEdge(netSol, netSol.Epmp);
    condEsp = (vertcat(Epmp.esp));
    
    dp = cell(numel(condEsp),1);
    if any(condEsp)
        Eesp = Epmp(condEsp);
        for i=1:numel(Eesp)
            vNet = getVertex(netSol, Eesp(i).vin);            
            vRes = getWellByName(wellSol,vNet.name);
            
            dp{i} = vNet.pressure-vRes.bhp;
        end
         try
             dp = vertcat(dp{:});
         catch
            error('Problem using function.  Assigning a value of 0.'); 
        end
    end
end



