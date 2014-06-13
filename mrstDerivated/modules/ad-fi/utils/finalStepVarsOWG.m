function obj = finalStepVarsOWG(step,finalState,wellSol,schedule,finalTime,fluid,system,varargin)
% the final state of the simulation (finalState) should equal stateNext

opt     = struct('ComputePartials',false,'xvScale',[],'xLeftSeed',[],'vLeftSeed',[]);

opt     = merge_options(opt, varargin{:});



dts   = schedule.step.val;
time = sum(dts(1:(step)));
if isfield(schedule,'time')
    time = time + schedule.time;
end

if finalTime ~= time
    
    % p has the appropriate dimension to instantiate a null vector
    p = finalState.pressure;
    qWs = vertcat(wellSol.qWs);

    
    if opt.ComputePartials
        [p, ~, ~, qWs, ~,~,~] = initVariablesADI(p,p,p,qWs,qWs,qWs,qWs); % place holders
    end
    
    p = 0*p;
    qWs = qWs*0;
    
    obj = [p; p; p; qWs; qWs; qWs,qWs];
else
    
    % The transformation function may be given as an input and
    % generalized    
    [ p,sW,rGH ] = stateMrst2statePsWrGH(finalState,fluid,system,'partials',opt.ComputePartials);
    
    
    qWs = vertcat(wellSol.qWs);
    qOs = vertcat(wellSol.qOs);
    qGs = vertcat(wellSol.qGs);
    pBH = vertcat(wellSol.bhp);
    
    
    if opt.ComputePartials
        % instantiating jacobians with right values and dimensions.
        % Observe that the independet variables are p,sw,x,qw,qo,qg,bhp
        % and not what the line below suggest!
        [pADI,sWADI,rGHADI,qWs,qOs,qGs,pBH] = initVariablesADI(double(p),double(sW),double(rGH),qWs,qOs,qGs,pBH);
        
        
        % establish jacobian given by stateMrst2statePsWrGH
        for k = 4:numel(pADI.jac)
            p.jac{k} = pADI.jac{k};
            sW.jac{k} = sWADI.jac{k};
            rGH.jac{k} = rGHADI.jac{k};
        end
        
    end
    
    obj = [p; sW; rGH; qWs;qOs; qGs; pBH];
    
    if ~isempty(opt.xvScale)
        obj = obj./[opt.xvScale];
    end
    
end



if ~isempty(opt.xLeftSeed)
    obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false);
end



end




