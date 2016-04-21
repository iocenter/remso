function netSol = runNetwork(ns, wellSol, forwardState,p, pScale,  varargin)
%% runNetwork: performs a network simulation of the production gathering network
    
    opt     = struct('dpFunction', @dpBeggsBrillJDJ, ...
                     'forwardGradient',true,...
                     'finiteDiff', false, ...
                     'ComputePartials',false,...
                     'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0), ...
                     'hasGas', false, ...
                     'fluid', []);                     
                 
    opt     = merge_options(opt, varargin{:});

    comp = opt.activeComponents;


if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)        
      ns = setWellSolValues(ns, wellSol, forwardState, p, pScale, 'ComputePartials',opt.ComputePartials, 'activeComponents', comp, 'hasGas', false, 'fluid', opt.fluid);        
elseif comp.gas && comp.oil && comp.water        
      ns = setWellSolValues(ns, wellSol, forwardState, p, pScale, 'ComputePartials',opt.ComputePartials, 'activeComponents', comp, 'hasGas', true, 'fluid', opt.fluid);                
else
    error('Not implemented for current activeComponents');
end

qoV = ns.qo(ns.VwProd);
qwV = ns.qw(ns.VwProd);
qgV = ns.qg(ns.VwProd);

% Propagates Flows in the Network
ns.qo = ns.M'*qoV;
ns.qw = ns.M'*qwV;
ns.qg = ns.M'*qgV;

%% source vertices
Vin = getVertex(ns, ns.Vsrc);

% sinks
surfaceSinks = setdiff(vertcat(ns.Vsnk), vertcat(ns.VwInj));
Vout = getVertex(ns, vertcat(surfaceSinks));

%% correction for multiple input pipelines
Ein = getEdge(ns, Vout.Ein);
Vout = getVertex(ns, Ein.vout);

%% propagate pressures in the network forward and backwards simultaneously
ns.pV = propagateDp(ns, Vin, Vout, opt.dpFunction, opt.hasGas, opt.forwardGradient, opt.finiteDiff);

netSol = ns;
end

function [pV] = propagateDp(ns, Vin, Vout, dpFunction, hasGas, forwardGradient, finiteDiff)
% propagateDp: calculates pressures of all nodes of the network.
Eout =  getEdge(ns, vertcat(Vin.Eout));   
Ein =  getEdge(ns, vertcat(Vout.Ein));  
condStopOut = vertcat(Eout.equipment);
condStopIn = vertcat(Ein.equipment);

while  ~all(condStopOut) || ~all(condStopIn)       
    Ef = [Eout(~condStopOut); Ein(~condStopIn)];    
    voutAux = getVertex(ns, vertcat(Ein(~condStopIn).vout));   %% correct size of output vertices  (manifold)    
    
    Vf = [Vin(~condStopOut)'; voutAux' ];
    negSign = ones(numel(Eout(~condStopOut)),1).*-1;
    posSign = ones(numel(Ein(~condStopIn)),1);
    sign = [negSign; posSign];
    
    [dp, pvK] = vectorizedDp(ns, Ef, Vf, dpFunction, hasGas, forwardGradient, finiteDiff);
  
    idsVout = vertcat(Eout(~condStopOut).vout);
    idsVin = vertcat(Ein(~condStopIn).vin);
    idsf = [idsVout; idsVin];    
    
    ns.pV(idsf) = pvK + sign.*dp;
    
    Vin = getVertex(ns, idsVout);
    Eout = getEdge(ns,  vertcat(Vin.Eout));
    
    if numel(Vin) ~= numel(Eout)
        error('Manifold before the equipment in the network or splitting of flows is happening');
    end
    
    Vout = getVertex(ns, idsVin);
    Ein = getEdge(ns,  vertcat(Vout.Ein));
    
    condStopOut = vertcat(Eout.equipment);     
    condStopIn = vertcat(Ein.equipment);             
end
pV = ns.pV;
end

function [pV] = forwardDp(ns, Vin, dpFunction, hasGas, forwardGradient, finiteDiff)
% forwardDp: performs foward pressure drop calculations in the network.
Eout =  getEdge(ns, vertcat(Vin.Eout));   
condStop = vertcat(Eout.equipment);
while  ~all(condStop)       
    [dp, pvK] = vectorizedDp(ns, Eout(~condStop), Vin(~condStop), dpFunction, hasGas, forwardGradient, finiteDiff);
  
    idsVout = vertcat(Eout(~condStop).vout);
    ns.pV(idsVout) = pvK - dp;
    
    Vin = getVertex(ns, idsVout);
    Eout = getEdge(ns,  vertcat(Vin.Eout));
    
    if numel(Vin) ~= numel(Eout)
        error('Manifold before the equipment in the network or splitting of flows is happening');
    end
    
    condStop = vertcat(Eout.equipment);     
end
pV = ns.pV;
end


function [pV] = backwardDp(ns, Vout, dpFunction, hasGas, forwardGradient, finiteDiff)
% forwardDp: performs backward pressure drop calculations in the network.
Ein =  getEdge(ns, vertcat(Vout.Ein));   
condStop = vertcat(Ein.equipment);

while  ~all(condStop)    
    voutAux = getVertex(ns, vertcat(Ein(~condStop).vout));   %% to correct size of output vertices  (manifold)    
    
    [dp, pvK] = vectorizedDp(ns, Ein(~condStop), voutAux, dpFunction,  hasGas, forwardGradient, finiteDiff);    
  
    idsVin = vertcat(Ein(~condStop).vin);    
    ns.pV(idsVin) = pvK + dp;
    
    Vout = getVertex(ns, idsVin);
    Ein = getEdge(ns,  vertcat(Vout.Ein));
    
    condStop = vertcat(Ein.equipment);            
end
pV = ns.pV;
end

function [dp, pvK] = vectorizedDp(ns, E, V, dpFunction, hasGas, forwardGradient, finiteDiff)
idsE = vertcat(E.id);        
idsV = vertcat(V.id);
qoK = ns.qo(idsE);
qwK = ns.qw(idsE);
qgK = ns.qg(idsE);

pvK = ns.pV(idsV);

dp = dpCVODES(E, qoK, qwK, qgK, pvK, 'dpFunction', dpFunction, 'hasSurfaceGas', hasGas,'forwardGradient', forwardGradient, 'finiteDiff', finiteDiff);
end



