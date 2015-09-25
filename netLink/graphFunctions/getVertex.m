function [ v ] = getVertex(ns, idV )
%GETVERTEX get the vertex mock object in the network object (ns) using its
%id (idV)
    V = ns.V;
    idx = ismember(vertcat(V.id), idV);
    if ~isempty(idx)
        v = V(idx);
    else
        v = [];
    end
end

