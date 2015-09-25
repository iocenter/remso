function [ e ] = getEdge(ns, idE)
%GETEDGE get the edge mock object in the network object (ns) using its id
%(idE)   
    E = ns.E;
    idx = ismember(vertcat(E.id), idE);
    if ~isempty(idx)
        e = E(idx);   
    else
        e = [];
    end
end

