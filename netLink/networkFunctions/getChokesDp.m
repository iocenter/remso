function [ dp] = getChokesDp(netSol)
%GETCHOKESDP get pressure drops in the network chokes
    Echk = getEdge(netSol, netSol.Eeqp);
    dp = cell(length(Echk),1);
    for i=1:length(Echk)
        vin = getVertex(netSol, Echk(i).vin);
        vout = getVertex(netSol, Echk(i).vout);
        
        dp{i} = vin.pressure-vout.pressure;
    end
    dp = vertcat(dp{:});
end

