function [ dp] = getChokesDp(netSol)
%GETCHOKESDP get pressure drops in the network chokes    
    Eeq = getEdge(netSol, netSol.Eeqp);
    dp = cell(length(Eeq),1);
    for i=1:length(Eeq)
        vin = getVertex(netSol, Eeq(i).vin);
        vout = getVertex(netSol, Eeq(i).vout);
        
        dp{i} = (vin.pressure-vout.pressure);
    end
    try
        dp = vertcat(dp{:});        
    catch
        error('Problem using function.  Assigning a value of 0.'); 
    end
end

