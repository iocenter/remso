function netSol = addEdge(ns, e, varargin)
% adds edge e to the network mock object ns
    opt     = struct('isEquipment', false); % default edge       
    opt     = merge_options(opt, varargin{:});
    
    % updating special edges, i.e., those leaving sources, sinks or
    % representing an equipement such as pumps, chokes or separators).
    if opt.isEquipment
         e.equipment = true;
         ns.Eeqp =  [ns.Eeqp; e.id];
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

