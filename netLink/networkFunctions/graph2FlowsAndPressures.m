function [ qo, qw, qg, p ] = graph2FlowsAndPressures(V,E)
%graph2FlowsAndPressures Creates flow and pressure vectors for a given
% graph G = (V,E).

    if length(V) ~= length(E)
        error(id('Graph:Pressure Drop Error'), ...
                   'Different sizes for edges and output vertices input to the pressure drop function');
    else                
        qo = abs(vertcat(E.qoE));
        qw = abs(vertcat(E.qwE));
        qg = abs(vertcat(E.qgE));
        p = vertcat(V.pressure);    
    end
end

