function objs = finalStepVarsOWG(forwardStates,schedule,fluid,activeComponents,varargin)
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


% The transformation function may be given as an input and
% generalized
disgas = activeComponents.disgas;
vapoil = activeComponents.vapoil;
[ p,sW,rGH ] = stateMrst2statePsWrGH(finalState,fluid,disgas,vapoil,'partials',opt.ComputePartials);
qWs = vertcat(wellSol.qWs);
qOs = vertcat(wellSol.qOs);
qGs = vertcat(wellSol.qGs);
pBH = vertcat(wellSol.bhp);


if opt.ComputePartials
    % instantiating jacobians with right values and dimensions.
    % Observe that the independet variables are p,sw,x,qw,qo,qg,bhp
    % and not what the line below suggest!
    [pADI,sWADI,rGHADI,qWs,qOs,qGs,pBH] = initVariablesADI(double(p),double(sW),double(rGH),qWs,qOs,qGs,pBH);
    
    
    % revert jacobian given by stateMrst2statePsWrGH
    for k = 4:numel(pADI.jac)
        p.jac{k} = pADI.jac{k};
        sW.jac{k} = sWADI.jac{k};
        rGH.jac{k} = rGHADI.jac{k};
    end
    
end


if K > 1
    
    p0 = p*0;
    q0 = pBH*0;
    
    obj0 = [p0;p0;p0;q0;q0;q0;q0];
    
    if opt.ComputePartials && ~(size(opt.xLeftSeed,2)==0)  %% to preserve the size
        obj0.jac = cellfun(@(x)sparse(size(opt.xLeftSeed,1),size(x,2)),obj0.jac,'UniformOutput',false);
    end
    
    for step = 1:K-1
        objs{step} = obj0;
    end
    
end

obj = [p; sW; rGH; qWs;qOs; qGs; pBH];

if ~isempty(opt.xvScale)
    obj = obj./[opt.xvScale];
end

if opt.ComputePartials && size(opt.xLeftSeed,2)>0
    obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false);
end

objs{K} = obj;


end






