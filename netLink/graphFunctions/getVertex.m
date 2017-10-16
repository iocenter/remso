function [ v ] = getVertex(ns, idV )
%GETVERTEX get the vertex mock object in the network object (ns) using its
%id (idV)
    v = ns.V(idV);
end

