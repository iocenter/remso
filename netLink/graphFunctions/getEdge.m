function [ e ] = getEdge(ns, idE)
%GETEDGE get the edge mock object in the network object (ns) using its id
%(idE)   
    e = ns.E(idE);    
end

