function [ M ] = createFlowMatrix( ns )
%createFlowMatrix Creates a flow matrix which contains a mapping of source
%nodes to edges of the gathering network

Vsrc = getVertex(ns, ns.Vsrc);
E = ns.E;
M = zeros(numel(Vsrc), numel(E));

for i=1:numel(Vsrc)
    Eout = getEdge(ns, Vsrc(i).Eout);
    
    while ~isempty(Eout)
        M(i,Eout.id) = 1;
        v = getVertex(ns, Eout.vout);
        
        Eout = getEdge(ns, v.Eout);
    end
end

end

