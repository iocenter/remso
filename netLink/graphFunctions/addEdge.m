function ns = addEdge(ns, e)
% adds edge e to the network mock object ns
    ns.E = [ns.E e];    
    ns.A(e.vin.id, e.vout.id) = e.id;
end

