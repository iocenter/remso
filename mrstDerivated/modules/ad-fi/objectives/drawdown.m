function obj = drawdown(forwardStates,schedule,fluid,varargin)

opt     = struct('ComputePartials',false,'leftSeed',[],'pscale',5*barsa);
opt     = merge_options(opt, varargin{:});


wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);


dts   = schedule.step.val;
numSteps = numel(dts);
tSteps = 1:numSteps;

obj = cell(1,numSteps);
for step = tSteps
    
    % pressure and saturaton vectors just used for place-holding
    p  = forwardStates{step}.pressure;
    sW = forwardStates{step}.s(:,1);
    wellSol = wellSols{step};
    qWs  = vertcat(wellSol.qWs);
    qOs  = vertcat(wellSol.qOs);
    pBHP = vertcat(wellSol.bhp);
    
    W = schedule.control(schedule.step.control(step)).W;    
    
    if opt.ComputePartials
        [p, sW, qWs, qOs, pBHP] = initVariablesADI(p, sW, qWs, qOs, pBHP);
    end
    
    g  = norm(gravity);
    [Tw, dzw, Rw, wc, perf2well, pInx, iInxW, iInxO] = getWellStuff(W);
    
    pw  = p(wc);
    
    
    bW     = fluid.bW(pw);
    rhoW   = bW.*fluid.rhoWS;
    
    bO     = fluid.bO(pw);
    rhoO   = bO.*fluid.rhoOS;
    
    
    dW  = (pBHP(perf2well) - pw + g*dzw.*rhoW)/opt.pscale;
    dO  = (pBHP(perf2well) - pw + g*dzw.*rhoO)/opt.pscale;
    
    d = pBHP;
    for w = 1:numel(W)
        perfw = find(Rw(:,w));
        d(w) = max(-W(w).sign*([dW(perfw);dO(perfw)]));
    end
    
    
    obj{step} = d;
    
    
    if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
        obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
    end
    
end
if numel(tSteps) > 1
    % set to zero intermediate values that are not the maximum
    [v,i] = max(cell2mat(cellfun(@(oi)double(oi),obj,'UniformOutput',false)),[],2);
    for step = tSteps
        W = schedule.control(schedule.step.control(step)).W;    
        for w = 1:numel(W)
            if i(w) ~= step
                obj{step}(w) = obj{step}(w)*0; 
            end
        end
    end
end