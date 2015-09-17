function ns = addVertex(ns, v)
% Adds vertex V to the network mock object ns
    ns.V = [ns.V v];
    B = zeros(size(ns.A,1)+1, size(ns.A,2)+1);
    for i=1:size(ns.A,1)
        for j=1:size(ns.A,2)
            B(i,j) = ns.A(i,j);
        end            
    end
    B(v.id, v.id) = v.id;
    ns.A = B;
end

