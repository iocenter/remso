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

% Forward propagation of pressures in the network
Vin = getVertex(ns, ns.Vsrc);
ns.pV = forwardDp(ns, Vin, 'dpFunction', opt.dpFunction, 'hasGas', opt.hasGas, 'forwardGradient', opt.forwardGradient, 'finiteDiff', opt.finiteDiff);

% Backward propagation of pressures in the network
surfaceSinks = setdiff(vertcat(ns.Vsnk), vertcat(ns.VwInj));
Vout = getVertex(ns, surfaceSinks);
ns.pV  = backwardDp(ns,Vout, 'dpFunction', opt.dpFunction, 'hasGas', opt.hasGas, 'forwardGradient', opt.forwardGradient, 'finiteDiff', opt.finiteDiff);

netSol = ns;
end

function [pV] = forwardDp(ns, Vin, varargin)
% forwardDp: performs foward pressure drop calculations in the network.
opt = struct('dpFunction', @dpBeggsBrillJDJ, ...
             'hasGas', false, ...
             'forwardGradient', true, ...
             'finiteDiff', true);
opt     = merge_options(opt, varargin{:});

Eout =  getEdge(ns, vertcat(Vin.Eout));   
condStop = vertcat(Eout.equipment);
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
    
    condStop = vertcat(Eout.equipment);    
end
pV = ns.pV;
end


function [pV] = backwardDp(ns, Vout, varargin)
% forwardDp: performs backward pressure drop calculations in the network.
opt = struct('dpFunction', @dpBeggsBrillJDJ, ...
             'hasGas', false, ...
             'forwardGradient', true, ...
             'finiteDiff', true);
opt     = merge_options(opt, varargin{:});

Ein =  getEdge(ns, vertcat(Vout.Ein));   
condStop = vertcat(Ein.equipment);

while  ~all(condStop)      
    [dp, pvK] = vectorizedBackwardDp(ns, Ein(~condStop), Vout(~condStop), 'dpFunction', opt.dpFunction, 'hasGas', opt.hasGas);    
  
    idsVin = vertcat(Ein(~condStop).vin);
    ns.pV(idsVin) = pvK + dp;
    
    Vout = getVertex(ns, idsVin);
    Ein = getEdge(ns,  vertcat(Vout.Ein));
    
    if numel(Vout) ~= numel(Ein) %% it is a manifold
        Vout = repmat(Vout, numel(Ein),1);
    end
    
    condStop = vertcat(Ein.equipment);            
end
pV = ns.pV;
end

function [dp, pvK] = vectorizedBackwardDp(ns, Ein, Vout, varargin)
opt = struct('dpFunction', @simpleDp, 'hasGas', false);
opt     = merge_options(opt, varargin{:});

idsEin = vertcat(Ein.id);    
idsVout  = vertcat(Vout.id);
    
qoK = ns.qo(idsEin);
qwK = ns.qw(idsEin);
qgK = ns.qg(idsEin);
pvK  = ns.pV(idsVout);

dp = dpCVODES(Ein, qoK, qwK, qgK, pvK, 'dpFunction', opt.dpFunction, 'hasSurfaceGas', opt.hasGas,'forwardGradient', true, 'finiteDiff', true);
end


