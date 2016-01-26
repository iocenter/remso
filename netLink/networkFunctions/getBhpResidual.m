function [ dp] = getBhpResidual(netSol)
%GETBHPRESIDUAL get bottom-hole pressure residual (difference from reservoir and network simulations)
    Eeqp = getEdge(netSol, netSol.Eeqp);        
    dp = cell(numel(Eeqp),1);        
    for i=1:numel(Eeqp)
        vRes = getVertex(netSol, Eeqp(i).vin);                        
        vNet = getVertex(netSol, Eeqp(i).vout);  

        dp{i} = vRes.pressure-vNet.pressure;
    end
     try
         dp = vertcat(dp{:});
     catch
        error('Problem using function.  Assigning a value of 0.'); 
    end
    
end



