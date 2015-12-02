function netSol = addEdge(ns, e, varargin)
% adds edge e to the network mock object ns
    opt     = struct('isEquipment', false, 'isChoke', false, 'isPump', false, 'isSeparator', false, 'isSource', false, 'isSink', false, 'isESP', false, 'isControllable', false); % default edge       
    opt     = merge_options(opt, varargin{:});
    
    % modifying adjacency matrix
    ns.A(e.vin, e.vout) = e.id;
    
    % updating special edges, i.e., those leaving sources, sinks or
    % representing an equipement such as pumps, chokes or separators).
    if opt.isEquipment % equipment edge (nodal analysis edge, for example)
         e.equipment = true;
         ns.Eeqp =  [ns.Eeqp; e.id];
    elseif opt.isChoke % production chokes
        e.equipment = true;
        e.choke = true;
        ns.Eeqp =  [ns.Eeqp; e.id];        
        ns.Echk =  [ns.Echk; e.id];        
    elseif opt.isPump  % pumps
        e.equipment = true;
        
        e.pump = opt.isPump;        
        
        ns.Eeqp =  [ns.Eeqp; e.id];        
        ns.Epmp =  [ns.Epmp; e.id];    
    elseif opt.isESP  % Special pumps for which an explicit model is known.
        e.equipment = false;
        e.esp  = opt.isESP;        
    elseif opt.isSeparator  % separators 
        e.equipment = true;
        e.separator = true;
        ns.Eeqp =  [ns.Eeqp; e.id];        
        ns.Esep =  [ns.Esep; e.id];    
    elseif opt.isSource % edges leaving a source node
        ns.Esrc = [ns.Esrc; e.id];        
    elseif opt.isSink % edges reaching a sink node
        ns.Esnk = [ns.Esnk; e.id];   
    elseif opt.isControllable
        ns.Ec = [ns.Ec; e.id];
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

