function grad = runGradientStep(G, rock, fluid, schedule, obj, S,W, varargin)


opt = struct('Verbose',    mrstVerbose, ...
    'ForwardStates',       [], ...
    'xRightSeeds',         1,...
    'uRightSeeds',         [],...
    'fwdJac',[]);

opt = merge_options(opt, varargin{:});


nAdj = size(obj.partials(2).s,1);
nFwd = size(opt.xRightSeeds,2);


assert(numel(opt.uRightSeeds)==1,'Not implemented for more');
assert(numel(opt.ForwardStates)==2)
assert(numel(obj.partials)==2)

if (nFwd > nAdj )
    controls = initControls(schedule);

    
    adjRes = runAdjoint(opt.ForwardStates, G, S, W, rock, fluid, schedule, controls,obj, 'Verbose',  opt.Verbose);
    
    
    [~,RHS] = solveAdjointTransportSystem(G, S, W, rock, fluid, opt.ForwardStates, adjRes, obj);
    
    gradx = -RHS;
    gradu = computeGradient(W, adjRes, schedule, controls);
    
    grad = gradu{1}'*opt.uRightSeeds{1}+gradx'*opt.xRightSeeds;
else
    
    
    controls = initControls(schedule);
    
    [grad ] = runForwardGradient(G, S, W, rock, fluid, opt.ForwardStates,schedule, controls,opt.xRightSeeds,opt.uRightSeeds{1}, obj);
    
end


end

