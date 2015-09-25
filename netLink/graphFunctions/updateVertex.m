function [ netSol ] = updateVertex(ns, vlist )
%UPDATEVERTEX update vertices in the network
    ids = vertcat(vlist.id);
    ns.V(ids) = vlist;
    netSol = ns;
end

