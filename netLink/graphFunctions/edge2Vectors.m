function [ qo, qw, qg, p] = edge2Vectors(E)
%EDGE2VECTORS creates flow and pressure vectors for a given set of edges 
% considering only output vertices and its pressures as decision variables.
    V = [];
    for i=1:length(E)
       V = [V; E(i).vout];
    end
    
    [qo, qw, qg, p] = graph2FlowsAndPressures(V,E);
end

