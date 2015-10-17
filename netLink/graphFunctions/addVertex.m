function ns = addVertex(ns, v, varargin)
% Adds vertex V to the network mock object ns
    opt     = struct('isProducer',false,'isInjector',false,'isSource',false,'isSink',false,...
        'isInterior', true, 'isControllable', false); % default vertex   
    
    opt     = merge_options(opt, varargin{:});

    ns.V = [ns.V v];
    
    % modifying adjacency matrix
    B = zeros(size(ns.A,1)+1, size(ns.A,2)+1);
    for i=1:size(ns.A,1)
        for j=1:size(ns.A,2)
            B(i,j) = ns.A(i,j);
        end            
    end
    B(v.id, v.id) = v.id;
    ns.A = B;
    
    % updating vertices sets
    if opt.isProducer
        ns.VwProd =  [ns.VwProd; v.id];  % producers
        ns.Vw     =  [ns.Vw; v.id];      % border node connecting with reservoir        
        ns.Vsrc   =  [ns.Vsrc; v.id];    % source node
    elseif opt.isInjector
        ns.VwInj  =  [ns.VwInj; v.id];   % injectors
        ns.Vw     =  [ns.Vw; v.id];      % bd node interfacing with reservoir
        ns.Vsnk   =  [ns.Vsnk; v.id];    % sink node
    elseif opt.isSource
        ns.Vsrc = [ns.Vsrc; v.id];       % source nodes
    elseif opt.isSink
        ns.Vsnk = [ns.Vsnk; v.id];       % sink nodes
    elseif opt.isInterior % interior node
        ns.Vint = [ns.Vint; v.id];       % interior nodes
    end
    
    if opt.isControllable
        ns.Vc = [ns.Vc; v.id];
    end
end

