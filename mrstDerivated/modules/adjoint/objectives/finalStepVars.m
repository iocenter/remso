function objs = finalStepVars(simRes, varargin)


opt     = struct('ComputePartials',false,'xvScale',[],'xLeftSeed',[],'vLeftSeed',[]);

opt     = merge_options(opt, varargin{:});


numSteps = numel(simRes);

objs = cell(1,numSteps);
step = numSteps;

finalState = simRes(step).resSol;

sW  = finalState.s(:,1);
q_w = vertcat(finalState.wellSol.flux);

if opt.ComputePartials
    [sW, q_w] = initVariablesADI(sW, q_w);
end

if numSteps > 1
    
    % set values and jacobians to zero
    sW0= 0*sW;
    q_W0 = 0*q_w;
    
    obj0 = [sW0; q_W0];
    
    if opt.ComputePartials && ~(size(opt.xLeftSeed,2)==0)  %% to preserve the size
        obj0.jac = cellfun(@(x)sparse(size(opt.xLeftSeed,1),size(x,2)),obj0.jac,'UniformOutput',false);
    end
    
    for step = 1:numSteps-1
        objs{step} = obj0;
    end
    
    
end

obj = [sW;q_w];

if ~isempty(opt.xvScale)
    obj = obj./[opt.xvScale];
end

if opt.ComputePartials && size(opt.xLeftSeed,2) > 0;
    obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false);
end

objs{numSteps} = obj;



end
