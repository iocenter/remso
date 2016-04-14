function netSol = runNetwork(ns, wellSol, forwardState,p, pScale,  varargin)
%% runNetwork: performs a network simulation of the production gathering network
    
    opt     = struct('dpFunction', @simpleDp, ...
                     'ComputePartials',false,...
                     'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0), ...
                     'hasGas', false, ...
                     'fluid', [],...
                     'sensitivityAnalysis', false, ...
                     'turnoffPumps', false);              
                 
    opt     = merge_options(opt, varargin{:});

    comp = opt.activeComponents;

if ~opt.sensitivityAnalysis
    if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)        
          ns = setWellSolValues(ns, wellSol, forwardState, p, pScale, 'ComputePartials',opt.ComputePartials, 'activeComponents', comp, 'hasGas', false, 'fluid', opt.fluid);        
    elseif comp.gas && comp.oil && comp.water        
          ns = setWellSolValues(ns, wellSol, forwardState, p, pScale, 'ComputePartials',opt.ComputePartials, 'activeComponents', comp, 'hasGas', true, 'fluid', opt.fluid);                
    else
        error('Not implemented for current activeComponents');
    end
end



%%TODO: qoV = flows of wells
qoV = ns.qo(ns.VwProd);
qwV = ns.qw(ns.VwProd);
qgV = ns.qg(ns.VwProd);

% Propagates Flows in the Network
ns.qo = ns.M'*qoV;
ns.qw = ns.M'*qwV;
ns.qg = ns.M'*qgV;

% Forward propagation of pressures in the network
Vin = getVertex(ns, ns.Vsrc);
ns.pV = forwardDp(ns, Vin,  'uptoChokeOrPump', true, 'dpFunction', opt.dpFunction, 'hasGas', opt.hasGas);

% Backward propagation of pressures in the network
surfaceSinks = setdiff(vertcat(ns.Vsnk), vertcat(ns.VwInj));
for j=1:numel(surfaceSinks)
    Vout = getVertex(ns, surfaceSinks(j));
    
    ns.pV  = backwardDp(ns,Vout, 'dpFunction', opt.dpFunction, 'hasGas', opt.hasGas);
end

netSol = ns;
end

function [pV] = forwardDp(ns, Vin, varargin)
% forwardDp: performs foward pressure drop calculations in the network.
opt = struct('uptoChokeOrPump', true, 'dpFunction', @simpleDp, 'hasGas', false);
opt     = merge_options(opt, varargin{:});

Eout =  getEdge(ns, vertcat(Vin.Eout));   
if opt.uptoChokeOrPump
    condStop = vertcat(Eout.equipment);
else
    condStop = ismember(vertcat(Eout.vout), vertcat(ns.Vsnk));
end
while  ~all(condStop)    
    % calculating pressure drops in the pipeline
    idsEout = vertcat(Eout.id);    
    idsVin  = vertcat(Vin.id);
    
    qoK = ns.qo(idsEout);
    qwK = ns.qw(idsEout);
    qgK = ns.qg(idsEout);
    pvK  = ns.pV(idsVin);

    dp = dpCVODES(Eout, qoK, qwK, qgK, pvK, 'dpFunction', opt.dpFunction, 'hasSurfaceGas', opt.hasGas,'forwardGradient', true, 'finiteDiff', true);      
    idsVout = vertcat(Eout.vout);
    ns.pV(idsVout) = pvK - dp;
    
    Vin = getVertex(ns, idsVout);
    Eout = getEdge(ns,  vertcat(Vin.Eout));
    
    if isempty(Eout)
        condStop = ones(length(Eout),1);
    elseif opt.uptoChokeOrPump
        condStop = vertcat(Eout.choke) | vertcat(Eout.pump);
    else
        condStop = zeros(numel(Eout),1);
    end
end
pV = ns.pV;
end


function [pV] = backwardDp(ns, Vout, varargin)
% forwardDp: performs backward pressure drop calculations in the network.
opt = struct('uptoChokeOrPump', true, 'dpFunction', @simpleDp, 'hasGas', false);
opt     = merge_options(opt, varargin{:});

Ein =  getEdge(ns, vertcat(Vout.Ein));   
if opt.uptoChokeOrPump
    condStop = vertcat(Ein.equipment);
else
    condStop = ismember(vertcat(Ein.vin), vertcat(ns.Vsrc));
end
while  ~all(condStop)    
    % calculating pressure drops in the pipeline
    idsEin = vertcat(Ein.id);    
    idsVout  = vertcat(Vout.id);
    
    qoK = ns.qo(idsEin);
    qwK = ns.qw(idsEin);
    qgK = ns.qg(idsEin);
    pvK  = ns.pV(idsVout);
    
    if numel(double(pvK)) ~= numel(double(qoK)) %% it is a manifold
        pvK = repmat(pvK,numel(double(qoK)),1);
    end        

    dp = dpCVODES(Ein, qoK, qwK, qgK, pvK, 'dpFunction', opt.dpFunction, 'hasSurfaceGas', opt.hasGas,'forwardGradient', true, 'finiteDiff', true);      
    idsVin = vertcat(Ein.vin);
    ns.pV(idsVin) = pvK + dp;
    
    Vout = getVertex(ns, idsVin);
    Ein = getEdge(ns,  vertcat(Vout.Ein));
    
    if isempty(Ein)
        condStop = ones(length(Ein),1);
    elseif opt.uptoChokeOrPump
        condStop = vertcat(Ein.choke) | vertcat(Ein.pump);
    else
        condStop = zeros(numel(Ein),1);
    end
end
pV = ns.pV;
end
