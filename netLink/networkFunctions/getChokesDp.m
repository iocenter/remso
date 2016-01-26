function [ dp] = getChokesDp(netSol)
%GETCHOKESDP get pressure drops in the network chokes
    if isempty(netSol.Echk)
        dp = 0;
        return;
    end
    Echk = getEdge(netSol, netSol.Echk); %TODO: consider both pumps and chokes
    dp = cell(length(Echk),1);
    for i=1:length(Echk)
        vin = getVertex(netSol, Echk(i).vin);
        vout = getVertex(netSol, Echk(i).vout);
        
        dp{i} = (vin.pressure-vout.pressure);
    end
    try
        dp = vertcat(dp{:});
    catch
        error('Problem using function.  Assigning a value of 0.'); 
    end
end

