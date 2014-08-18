function objs = finalStepVarsOW(forwardStates,schedule,finalTime, varargin)
% the final state of the simulation (finalState) should equal stateNext

opt     = struct('ComputePartials',false,'xvScale',[],'xLeftSeed',[],'vLeftSeed',[]);

opt     = merge_options(opt, varargin{:});


K = numel(forwardStates);
objs = cell(1,K);

dts   = schedule.step.val;
dtFrac = schedule.step.val/(sum(schedule.step.val));

for step = 1:K
    
    finalState = forwardStates{step};
    wellSol = forwardStates{step}.wellSol;
    
    p   = finalState.pressure;
    sW  = finalState.s(:,1);
    qWs = vertcat(wellSol.qWs);
    qOs = vertcat(wellSol.qOs);
    pBH = vertcat(wellSol.bhp);
    
    
    if opt.ComputePartials
        [p, sW, qWs, qOs, pBH] = initVariablesADI(p, sW, qWs, qOs, pBH);
    end
    
    
    time = sum(dts(1:(step)));
    if isfield(schedule,'time')
        time = time + schedule.time;
    end
    
    if finalTime ~= time
        p = 0*p;
        sW = 0*sW;
    end
    
    pBH = dtFrac(step)*pBH;
    qWs = dtFrac(step)*qWs;
    qOs = dtFrac(step)*qOs;
    
    
    obj = [p; sW; qWs; qOs; pBH];
    
    if ~isempty(opt.xvScale)
        obj = obj./[opt.xvScale];
    end
    
    if ~isempty(opt.xLeftSeed)
        obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false);
    end
    
    objs{step} = obj;
    
end

end




