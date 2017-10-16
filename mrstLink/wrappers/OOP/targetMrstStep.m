function [ varargout ] = targetMrstStep(x0,u,target,simulator,schedule,reservoirP,varargin)
%
% Simulate a single step and apply a target function on the results
%
% SYNOPSIS:
%  [f,Jac,convergence,simVars] = targetMrstStep(x0,u,target,simulator,schedule,reservoirP)
%  [f,Jac,convergence,simVars] = targetMrstStep(x0,u,target,simulator,schedule,reservoirP, 'pn', pv, ...)
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
%   (uScale) - control variables scaling
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
opt = struct('gradients',false,'uScale',[],'xRightSeeds',[],'uRightSeeds',[],'guessX',[],'guessV',[],'saveJacobians',true,'simVars',[], 'fixedWells', []);
opt = merge_options(opt, varargin{:});

nx = numel(x0);
nu = numel(u);
model = reservoirP.model;

simulate = true;  % if the simVars provided, skip the simulation part
if ~isempty(opt.simVars)
   % simulate = ~opt.simVars.convergence.converged;
   simulate = false;
end


if opt.gradients
    
    if (size(opt.uRightSeeds,1)==0)
        uRightSeeds = [speye(nu),sparse(nu,nx)];
        xRightSeeds = [sparse(nx,nu),speye(nx)];
    else
        uRightSeeds = opt.uRightSeeds;
        xRightSeeds = opt.xRightSeeds;
    end  
else
    uRightSeeds = [];
end


if simulate
    if opt.gradients
        [ shootingVars.state0,JacTX] = model.toMRSTStates(x0);
    else
        [ shootingVars.state0] = model.toMRSTStates(x0);   
    end
    [ shootingVars.schedule,uRightSeeds ] = controls2Schedule( u,schedule,'uScale',opt.uScale,...
    'partials',opt.gradients,'uRightSeeds',uRightSeeds, 'fixedWells', opt.fixedWells);
    
    % The guess is only given for the last simulation step.  Do something if there is any intermediate. 
    if ~isempty(opt.guessX)
        
        nScheduleSteps = numel(shootingVars.schedule.step.val);
        shootingGuess = cell(nScheduleSteps,1);
        
        [ shootingGuess{nScheduleSteps} ] = model.toMRSTStates( opt.guessX);
        
        if ~isempty(opt.guessV)
            
            W = shootingVars.schedule.control(end).W;
            [shootingGuess{nScheduleSteps}.wellSol] = model.toWellSol( opt.guessV,W,shootingGuess{nScheduleSteps});
            
        end
    else
        shootingGuess = [];
    end
    
    
    [shootingSol,JacRes,convergence] = simulator(shootingVars,reservoirP,'shootingGuess',shootingGuess);
    
    
    forwardStates = shootingSol.ForwardStates;
    scheduleSol = shootingSol.schedule;
    
    targetObjs = callArroba(target,{forwardStates,...
        scheduleSol},'ComputePartials', opt.gradients);
    
    simVars.forwardStates = forwardStates;
    simVars.schedule = shootingSol.schedule;
    if opt.saveJacobians
        simVars.JacRes = JacRes;
    else
        simVars.JacRes = [];
    end
    simVars.convergence = convergence;
	if opt.saveJacobians
        simVars.targetObjs = targetObjs;
	else
        simVars.targetObjs = cellfun(@(obj)double(obj),targetObjs,'UniformOutput',false);
	end
    
    
else
    
    forwardStates = opt.simVars.forwardStates;
    scheduleSol = opt.simVars.schedule;
    JacRes = opt.simVars.JacRes;   %% does this really worth?
    convergence = opt.simVars.convergence;
    targetObjs = opt.simVars.targetObjs;
    simVars = opt.simVars;
end



targetK = cellfun(@(obj)double(obj),targetObjs,'UniformOutput',false);
sumTarget = sum(cat(3,targetK{:}),3);

Jac = [];
if opt.gradients
    
    if ~simulate
        [ ~,uRightSeeds ] = controls2Schedule( u,schedule,...
            'uScale',opt.uScale,...
            'partials',opt.gradients,'uRightSeeds',uRightSeeds, 'fixedWells', opt.fixedWells);
        [ shootingVars.state0,JacTX] = model.toMRSTStates(x0);
   
        targetObjs = callArroba(target,{forwardStates,...
            scheduleSol},'ComputePartials', opt.gradients);
    end
    
    
    % unpack and group the left jacobians;
    lS = cellfun(@(x)cat(x),targetObjs,'UniformOutput',false);
    lS = cellfun(@(x)x.jac{1},lS,'UniformOutput',false);
    lSF =@(step) lS{step};
    
    
    % TODO: Include in Jac*Vector in transformation functions
    % there should be some little advantage!
    xRightSeeds = JacTX*xRightSeeds;
    
    gradients = computeGradientAD(shootingVars.state0, forwardStates, model, schedule, lSF,xRightSeeds,uRightSeeds);
    
    
    if (size(opt.uRightSeeds,1)==0)
        Jac.Ju = gradients(:,1:nu);
        Jac.Jx = gradients(:,nu+1:nu+nx);
    else
        Jac.J = gradients;
    end
   
end

varargout{1} = sumTarget;
varargout{2} = Jac;
varargout{3} = convergence;
varargout{4} = simVars;


end