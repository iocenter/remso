function [ varargout ] = targetMrstStep(x0,u,target,simulator,wellSol,schedule,reservoirP,varargin)
%
% Simulate a single step and apply a target function on the results
%
% SYNOPSIS:
%  [f,Jac,convergence,simVars] = targetMrstStep(x0,u,target,simulator,wellSol,schedule,reservoirP)
%  [f,Jac,convergence,simVars] = targetMrstStep(x0,u,target,simulator,wellSol,schedule,reservoirP, 'pn', pv, ...)
%
% PARAMETERS:
%
%   x0 - initial state in Remso format
%
%   u - controls in Remso format
%
%   target -  Function follwing structure defined in dummyMrstFunc
%
%   simulator - simulator function
%
%   wellSol - wellSol mock object
%
%   schedule - schedule mock object
%
%   reservoirP - reservoir parameters
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   gradients - true if partial derivatives are computed.
%
%   leftSeed - vector for vector-Jacobian product.
%
%   (xRightSeeds,uRightSeeds) - for Jacobian-vector product
%
%   (xScale,vScale,uScale) - Variables scaling
%
%   (guessX,guessV) - guess of the simulation result
%
%    simVars - for hot starting
%
% RETURNS:
%
%   f - value of the target function.
%
%   Jac - Jacobain of the targe function
%
%   convergence - convergence information of the simulation step
%
%   simVars - simulation variables
%
% SEE ALSO:
%
%
opt = struct('gradients',false,'xScale',[],'vScale',[],'uScale',[],'xRightSeeds',[],'uRightSeeds',[],'guessX',[],'guessV',[],'saveJacobians',true,'simVars',[]);
opt = merge_options(opt, varargin{:});

nx = numel(x0);
nu = numel(u);
nw = numel(schedule.control(1).W);
np = nu-nw;


simulate = true;  % if the simVars provided, skip the simulation part
if ~isempty(opt.simVars)
   % simulate = ~opt.simVars.convergence.converged;
   simulate = false;
end


if opt.gradients
    
    if (size(opt.uRightSeeds,1)==0)
        wRightSeeds = [speye(nw),sparse(nw,np),sparse(nw,nx)];
        pRightSeeds = [sparse(np,nw),speye(np),sparse(np,nx)];
        xRightSeeds = [sparse(nx,nw),sparse(nx,np),speye(nx)];
    else
        wRightSeeds = opt.uRightSeeds(1:nw,:);
        pRightSeeds = opt.uRightSeeds(nw+1:end,:);
        xRightSeeds = opt.xRightSeeds;
    end  
else
    wRightSeeds = [];
end
w = u(1:nw);
p = u(nw+1:end);



if simulate
    
    [ shootingVars.state0,JacTX] = stateVector2stateMrst( x0,'xScale',opt.xScale,...
        'activeComponents',reservoirP.system.activeComponents,...
        'fluid',reservoirP.fluid,...
        'partials',opt.gradients);
    [ shootingVars.schedule,wRightSeeds ] = controls2Schedule( w,schedule,'uScale',opt.uScale,...
    'partials',opt.gradients,'uRightSeeds',wRightSeeds);
    
    % The guess is only given for the last simulation step.  Do something if there is any intermediate. 
    if ~isempty(opt.guessX)
        
        nScheduleSteps = numel(shootingVars.schedule.step.val);
        shootingGuess = cell(nScheduleSteps,1);
        [ shootingGuess{nScheduleSteps} ] = stateVector2stateMrst( opt.guessX,'xScale',opt.xScale,...
            'activeComponents',reservoirP.system.activeComponents,...
            'fluid',reservoirP.fluid);
        if ~isempty(opt.guessV)
            [shootingGuess{nScheduleSteps}.wellSol] = algVar2wellSol( opt.guessV,wellSol,'vScale',opt.vScale,...
                'activeComponents',reservoirP.system.activeComponents);
        end
    else
        shootingGuess = [];
    end
    
    
    [shootingSol,JacRes,convergence] = callArroba(simulator,{shootingVars,reservoirP},'shootingGuess',shootingGuess);
    
    
    forwardStates = shootingSol.ForwardStates;
    scheduleSol = shootingSol.schedule;
    
    targetObjs = callArroba(target,{forwardStates,...
        scheduleSol,p},'ComputePartials', opt.gradients);
    
    simVars.forwardStates = forwardStates;
    simVars.schedule = shootingSol.schedule;
    if opt.saveJacobians
        simVars.JacRes = JacRes;
    else
        simVars.JacRes = [];
    end
    simVars.convergence = convergence;
    simVars.targetObjs = cellfun(@(obj)double(obj),targetObjs,'UniformOutput',false);
	
  
else
    
    forwardStates = opt.simVars.forwardStates;
    scheduleSol = opt.simVars.schedule;
    JacRes = opt.simVars.JacRes;
    convergence = opt.simVars.convergence;
    targetObjs = opt.simVars.targetObjs;
    simVars = opt.simVars;
end



targetK = cellfun(@(obj)double(obj),targetObjs,'UniformOutput',false);
sumTarget = sum(cat(3,targetK{:}),3);

Jac = [];
if opt.gradients
    
    if ~simulate
        [ ~,wRightSeeds ] = controls2Schedule( w,schedule,...
            'uScale',opt.uScale,...
            'partials',opt.gradients,'uRightSeeds',wRightSeeds);
        [ shootingVars.state0,JacTX ] = stateVector2stateMrst( x0,...
            'xScale',opt.xScale,...
            'activeComponents',reservoirP.system.activeComponents,...
            'fluid',reservoirP.fluid,...
            'partials',opt.gradients);
        targetObjs = callArroba(target,{forwardStates,...
            scheduleSol,p},'ComputePartials', opt.gradients);
    end        
    Jacp = targetObjs{1}.jac{end};
    targetObjs{1} = ADI(targetObjs{1}.val,targetObjs{1}.jac(1:end-1));
    for j = 2:numel(targetObjs)
        Jacp = Jacp + targetObjs{j}.jac{end};
        targetObjs{j} = ADI(targetObjs{j}.val,targetObjs{j}.jac(1:end-1));
    end    

%     Jacp = targetObjs{1}.jac{6};
%     targetObjs{1} = ADI(targetObjs{1}.val,targetObjs{1}.jac(1:5));
%     for j = 2:numel(targetObjs)
%         Jacp = Jacp + targetObjs{j}.jac{6};
%         targetObjs{j} = ADI(targetObjs{j}.val,targetObjs{j}.jac(1:5));
%     end
    
    
    % unpack and group the left jacobians;
    lS = cellfun(@(x)cat(x),targetObjs,'UniformOutput',false);
    lS = cellfun(@(x)x.jac{1},lS,'UniformOutput',false);
    lSF =@(step) lS{step};
    
    
    % TODO: Include in Jac*Vector in transformation functions
    % there should be some little advantage!
    xRightSeeds = JacTX*xRightSeeds;
    
    gradients = runGradientStep(reservoirP.G, ...
        reservoirP.rock, ...
        reservoirP.fluid, ...
        scheduleSol,...
        lSF,...
        reservoirP.system,...
        'Verbose', false,...
        'ForwardStates', [{shootingVars.state0},forwardStates],...
        'xRightSeeds',xRightSeeds,...
        'uRightSeeds',wRightSeeds,...
        'fwdJac',JacRes);
    
    
    if (size(opt.uRightSeeds,1)==0)
        Jac.Ju = [gradients(:,1:nw),Jacp];
        Jac.Jx = gradients(:,nw+np+1:end);
    else
        Jac.J = gradients+Jacp*pRightSeeds;
    end
   
end

varargout{1} = sumTarget;
varargout{2} = Jac;
varargout{3} = convergence;
varargout{4} = simVars;


end