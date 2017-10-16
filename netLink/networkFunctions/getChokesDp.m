function [ dp] = getChokesDp(netSol)
%GETCHOKESDP get pressure drops in the network chokes    
    if isempty(netSol.Eeqp)
       error('There is no equipment in the network to obtain the pressure difference' )
    end
    Eeq = getEdge(netSol, netSol.Eeqp);
    dp = cell(length(Eeq),1);
    for i=1:length(Eeq)        
        dp{i} = netSol.pV(Eeq(i).vin)-netSol.pV(Eeq(i).vout); 
    end
    try
        dp = vertcat(dp{:});        
    catch
        error('Problem using function.  Assigning a value of 0.'); 
    end
end

