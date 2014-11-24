function [x,v,Jac,convergence,simVars] = mrstStep(x0,u,simulator,wellSol,schedule,reservoirP,varargin)
%
% Simulate a single step and apply a target function on the results
%
% SYNOPSIS:
%  [x,v,Jac,convergence,simVars] = mrstStep(x0,u,simulator,wellSol,schedule,reservoirP)
%  [x,v,Jac,convergence,simVars] = mrstStep(x0,u,simulator,wellSol,schedule,reservoirP, 'pn', pv, ...)
%
% PARAMETERS:
%
%   x0 - initial state in Remso format
%
%   u - controls in Remso format
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
%   (xLeftSeed,vLeftSeed) - vector for vector-Jacobian product.
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
%   xs - states at the end of the step.
%
%   vs - algebraic states at the end of the step
%
%   Jac - Jacobian
%
%   convergence - convergence information of the simulation step
%
%   simVars - simulation variables
%
% SEE ALSO:
%
%
opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'guessX',[],'guessV',[],'xScale',[],'vScale',[],'uScale',[],'saveJacobians',true,'simVars',[]);
opt = merge_options(opt, varargin{:});

finalTime = sum(schedule.step.val);
if isfield(schedule,'time')
    finalTime = finalTime + schedule.time;
end
   
nx = numel(opt.xScale);
nv = numel(opt.vScale);
nw = numel(wellSol);

target =@(j,shootingSolN,wellSol,schedule,varargin) finalStepVars(j,shootingSolN,wellSol,schedule,finalTime,'xvScale',[opt.xScale;opt.vScale],...
                                                                  'xLeftSeed',opt.xLeftSeed,'vLeftSeed',opt.vLeftSeed,varargin{:});
% TODO: Change for OWG
if ~isempty(opt.guessV)
   opt.guessV = opt.guessV(1:3*nw);
end
                                                              
[f,J,convergence,simVars] = targetMrstStep(x0,u,target,simulator,wellSol,schedule,reservoirP,...
    'gradients',opt.gradients,...
    'xScale',opt.xScale,...
    'vScale',opt.vScale(1:3*nw),...
    'uScale',opt.uScale,...
    'xRightSeeds',opt.xRightSeeds,...
    'uRightSeeds',opt.uRightSeeds,...
    'guessX',opt.guessX,...
    'guessV',opt.guessV,...
    'saveJacobians',opt.saveJacobians,...
    'simVars',opt.simVars);

x = f(1:nx);
v = f(nx+1:nx+nv);

Jac = [];
if opt.gradients 
    if ~(size(opt.xLeftSeed,2)==0) && ~(size(opt.xRightSeeds,1)==0)
        error('not implemented')
    elseif (size(opt.xLeftSeed,2)==0) && (size(opt.xRightSeeds,1)==0)
        Jac.xJu  = J.Ju(1:nx,:);
        Jac.vJu  = J.Ju(nx+1:nx+nv,:);
        Jac.xJx  = J.Jx(1:nx,:);
        Jac.vJx  = J.Jx(nx+1:nx+nv,:);
    elseif ~(size(opt.xLeftSeed,2)==0)
        Jac.Ju  = J.Ju;
        Jac.Jx  = J.Jx;
	elseif ~(size(opt.xRightSeeds,1)==0)
        Jac.xJ  = J.J(1:nx,:);
        Jac.vJ  = J.J(nx+1:nx+nv,:);
    else
        error('what')
    end
end



end
