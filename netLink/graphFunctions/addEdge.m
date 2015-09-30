function netSol = addEdge(ns, e, varargin)
% adds edge e to the network mock object ns
    opt     = struct('isEquipment',false, 'isSource', false, 'isSink', false); % default edge       
    opt     = merge_options(opt, varargin{:});

   
    
    % modifying adjacency matrix
    ns.A(e.vin, e.vout) = e.id;
    
    % updating edges sets   
    if opt.isEquipment   % special equipment such chokes, compressors
        e.equipment = true;
        ns.Eeqp =  [ns.Eeqp; e.id];        
    elseif opt.isSource % edges leaving a source node
        ns.Esrc = [ns.Esrc; e.id];        
    elseif opt.isSink % edges reaching a sink node
        ns.Esnk = [ns.Esnk; e.id];
    end
    
    % adding new edge to the main edges list
    ns.E = [ns.E; e];  
    
    % updating vertices affected with the addition of the edge 
    vorig = getVertex(ns, e.vin);
    vorig.Eout = [vorig.Eout; e.id];
    
    vdest  = getVertex(ns, e.vout);
    vdest.Ein = [vdest.Ein; e.id];
    
    ns = updateVertex(ns, [vorig; vdest]);    
    
    netSol = ns;
end

