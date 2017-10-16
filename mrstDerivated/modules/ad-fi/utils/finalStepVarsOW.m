function objs = finalStepVarsOW(forwardStates,schedule,p, varargin)
% Final state of the simulation togehter with the algebraic variables at
% the end of the integration period.

%{
This function is designed to work togehter with targetMrstStep and
runGradientStep.  Therefore the size of objs depend on the number of steps
within the simulation schedule.

According to the implementation of runGradientStep, the output is the sum
of all objs, therefore only the last values must remain non-zero.

%}

opt     = struct('ComputePartials',false,'xvScale',[],'xLeftSeed',[],'vLeftSeed',[]);

opt     = merge_options(opt, varargin{:});


K = numel(forwardStates);
objs = cell(1,K);

%dts   = schedule.step.val;
%dtFrac = schedule.step.val/(sum(schedule.step.val));


step = K;

finalState = forwardStates{step};
wellSol = forwardStates{step}.wellSol;

pressure   = finalState.pressure;
sW  = finalState.s(:,1);
qWs = vertcat(wellSol.qWs);
qOs = vertcat(wellSol.qOs);
pBH = vertcat(wellSol.bhp);


if opt.ComputePartials
    [pressure, sW, qWs, qOs, pBH,~] = initVariablesADI(pressure, sW, qWs, qOs, pBH,p);
end


if K > 1
    
    % set values and jacobians to zero
    pressure0 = 0*pressure;
    sW0= 0*sW;
    pBH0 = 0*pBH;
    qWs0 = 0*qWs;
    qOs0 = 0*qOs;
    
    obj0 = [pressure0; sW0; qWs0; qOs0; pBH0];
    
    if opt.ComputePartials && ~(size(opt.xLeftSeed,2)==0)  %% to preserve the size
        obj0.jac = cellfun(@(x)sparse(size(opt.xLeftSeed,1),size(x,2)),obj0.jac,'UniformOutput',false);
    end
    
    for step = 1:K-1
        objs{step} = obj0;
    end
    
    
end

obj = [pressure; sW; qWs; qOs; pBH];

if ~isempty(opt.xvScale)
    obj = obj./[opt.xvScale];
end

if opt.ComputePartials && size(opt.xLeftSeed,2) > 0;
    obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false);
end

objs{K} = obj;


end




