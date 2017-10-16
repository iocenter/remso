function ns = addVertex(ns, v, varargin)
% Adds vertex V to the network mock object ns
    opt     = struct('isProducer',false,'isInjector',false,'isSource',false,'isSink',false); % default vertex   
    opt     = merge_options(opt, varargin{:});

    ns.V = [ns.V v];
    
    % updating vertices sets
    if opt.isProducer
        ns.VwProd =  [ns.VwProd; v.id];  % producers
        ns.Vw     =  [ns.Vw; v.id];      % border node connecting with reservoir        
        ns.Vsrc   =  [ns.Vsrc; v.id];    % source node
    elseif opt.isInjector
        ns.VwInj  =  [ns.VwInj; v.id];   % injectors
        ns.Vw     =  [ns.Vw; v.id];      % border nodes in the interface of network and reservoir
        ns.Vsnk   =  [ns.Vsnk; v.id];    % sink node
    end   
    
    if opt.isSource
        ns.Vsrc = [ns.Vsrc; v.id];       % source nodes
    elseif opt.isSink
        ns.Vsnk = [ns.Vsnk; v.id];       % sink nodes
    end

end

