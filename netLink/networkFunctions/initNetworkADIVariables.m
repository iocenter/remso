function [netSol] = initNetworkADIVariables(netSol, varargin )
%INITNETWORKADIVARIABLES Initializes Network ADI Variables
    E = vertcat(netSol.E);
    V = vertcat(netSol.V);

    [E.qoE, E.qwE, E.qgE, V.pressure] = initVariablesADI(E.qoE, E.qwE, E.qgE, V.pressure);
    
    netSol.E = E;
    netSol.V = V;
end

