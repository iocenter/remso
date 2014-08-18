function objs = finalStepVarsOWG(forwardStates,schedule,finalTime,fluid,system,varargin)
% the final state of the simulation (finalState) should equal stateNext
% wellsols is the weighted sum in time

opt     = struct('ComputePartials',false,'xvScale',[],'xLeftSeed',[],'vLeftSeed',[]);

opt     = merge_options(opt, varargin{:});


Ks = numel(forwardStates);

objs = cell(1,Ks);

dts   = schedule.step.val;
dtFrac = schedule.step.val/(sum(schedule.step.val));

for step = 1:Ks
    
    finalState = forwardStates{step};
    wellSol = forwardStates{step}.wellSol;
    
    
    time = sum(dts(1:(step)));
    if isfield(schedule,'time')
        time = time + schedule.time;
    end
    
    if finalTime ~= time
        
        % p has the appropriate dimension to instantiate a null vector
        p = finalState.pressure;
        
        
        qWs = vertcat(wellSol.qWs);
        qOs = vertcat(wellSol.qOs);
        qGs = vertcat(wellSol.qGs);
        pBH = vertcat(wellSol.bhp);
        
        
        if opt.ComputePartials
            [p, ~, ~,qWs,qOs,qGs,pBH] = initVariablesADI(p,p,p,qWs,qOs,qGs,pBH); % place holders
        end
        
        qWs = qWs*dtFrac(step);
        qOs = qOs*dtFrac(step);
        qGs = qGs*dtFrac(step);
        pBH = pBH*dtFrac(step);
        
        p = 0*p;
        
        obj = [p; p; p;qWs;qOs;qGs;pBH];
        
        if ~isempty(opt.xvScale)
            obj = obj./[opt.xvScale];
        end
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
        
        qWs = qWs*dtFrac(step);
        qOs = qOs*dtFrac(step);
        qGs = qGs*dtFrac(step);
        pBH = pBH*dtFrac(step);
        
        obj = [p; sW; rGH; qWs;qOs; qGs; pBH];
        
        if ~isempty(opt.xvScale)
            obj = obj./[opt.xvScale];
        end
        
    end
    
    if ~isempty(opt.xLeftSeed)
        obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false);
    end
    
    objs{step} = obj;
    
end


end




