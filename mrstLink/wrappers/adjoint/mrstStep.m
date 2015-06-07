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
opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'guessX',[],'guessV',[],'xScale',[],'vScale',[],'uScale',[],'saveJacobians',true,'simVars',[],'algFun',[]);
opt = merge_options(opt, varargin{:});

nx = reservoirP.G.cells.num;

nw = numel(wellSol);
nvw = sum(arrayfun(@(wi)numel(wi.cells),reservoirP.W));


if opt.gradients && ~(size(opt.vLeftSeed,2)==0)
    vwLeftSeed = opt.vLeftSeed(:,1:nvw);
    vnLeftSeed = opt.vLeftSeed(:,nvw+1:end);
    sumLeftSeeds = true;
else
    vwLeftSeed = [];
    vnLeftSeed = [];
    sumLeftSeeds = false;
end

target1 = arroba(@finalStepVars,[1,-1],{'xvScale',[opt.xScale;opt.vScale(1:nvw)],...
    'xLeftSeed',opt.xLeftSeed,'vLeftSeed',vwLeftSeed},...
    true);



if ~isempty(opt.algFun) %% merge the targets
    
    target2 = arroba(opt.algFun,[1,2],{'leftSeed',vnLeftSeed},...
        true);
    
    
    [ target ] = concatenateMrstTargets([target1,target2],sumLeftSeeds);
else
    
    target = target1;
    
end

[f,J,convergence,simVars] = targetMrstStep(x0,u,target,simulator,wellSol,schedule,reservoirP,...
    'gradients',opt.gradients,...
    'xScale',opt.xScale,...
    'vScale',opt.vScale,...
    'uScale',opt.uScale,...
    'xRightSeeds',opt.xRightSeeds,...
    'uRightSeeds',opt.uRightSeeds,...
    'guessX',opt.guessX,...
    'guessV',opt.guessV,...
    'saveJacobians',opt.saveJacobians,...
    'simVars',opt.simVars);

x = f(1:nx);
v = f(nx+1:end);

Jac = [];
if opt.gradients
    if ~(size(opt.xLeftSeed,2)==0) && ~(size(opt.xRightSeeds,1)==0)
        error('not implemented')
    elseif (size(opt.xLeftSeed,2)==0) && (size(opt.xRightSeeds,1)==0)
        Jac.xJu  = J.Ju(1:nx,:);
        Jac.vJu  = J.Ju(nx+1:end,:);
        Jac.xJx  = J.Jx(1:nx,:);
        Jac.vJx  = J.Jx(nx+1:end,:);
    elseif ~(size(opt.xLeftSeed,2)==0)
        Jac.Ju  = J.Ju;
        Jac.Jx  = J.Jx;
    elseif ~(size(opt.xRightSeeds,1)==0)
        Jac.xJ  = J.J(1:nx,:);
        Jac.vJ  = J.J(nx+1:end,:);
    else
        error('what')
    end
end



end
