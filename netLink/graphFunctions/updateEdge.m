function [ netSol ] = updateEdge(ns, elist)
%UPDATEEDGE Update edges in the network
    ids = vertcat(elist.id);
    ns.E(ids) = elist;
    netSol = ns;
end

