function plotDownstreamFlow(netSol, idV)    
    vi = getVertex(netSol, idV);
    while ~isempty(vi.Eout)
        fprintf('%s [%d] : %d \n', vi.name, vi.id, vi.pressure);
        ei = getEdge(netSol, vi.Eout);
        vi = getVertex(netSol, ei.vout);
    end
    fprintf('%s [%d] : %d \n', vi.name, vi.id, vi.pressure);
end