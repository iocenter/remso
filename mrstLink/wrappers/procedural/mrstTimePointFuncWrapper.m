function [f,Jac] = mrstTimePointFuncWrapper(xfk,uk,vk,target,schedule,wellSol,netSol,fluid,system,varargin)
%
% Interface with Remso a mrstPointFunction
%
% SYNOPSIS:
%  [f,Jac] = mrstTimePointFuncWrapper(xfk,uk,vk,target,schedule,wellSol,activeComponents,)
%  [f,Jac] = mrstTimePointFuncWrapper(xfk,uk,vk,target,schedule,wellSol,activeComponents, 'pn', pv, ...)
%
% PARAMETERS:
%
%   xfk - state
%
%   uk - control
%
%   vfk -  algebraic state.
%
%   target -  Function follwing structure defined in dummyMrstFunc
%
%   schedule - schedule mock object
%
%   wellSol - wellSol mock object
%
%   fluid      - Fluid as defined by initDeckADIFluid.
%
%   system     - System configuration as defined by initADISystem.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   partials - true if partial derivatives are computed.
%
%   leftSeed - vector for vector-Jacobian product.
%
%   (xRightSeeds,uRightSeeds,vRightSeeds) - for Jacobian-vector product
%
%   (xScale,vScale,uScale) - Variable scaling
%
% RETURNS:
%
%   f - value of the target function.
%
%   Jac - Jacobain of the target function
%
% SEE ALSO:
%
%

%
%TODO:  Deprecate, talk with Thiago

opt = struct('partials',false,'leftSeed',[],'xScale',[],'vScale',[],'uScale',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});

[wellSol,netSol,JacTW,JacTN] = algVar2mrstAlg(vk,wellSol,netSol,'vScale',opt.vScale,'partials',opt.partials,'activeComponents',system.activeComponents);
[ state,JacTX] = stateVector2stateMrst( xfk,'xScale',opt.xScale,...
    'activeComponents',system.activeComponents,...
    'fluid',fluid,...
    'system',system,...
    'partials',opt.partials);
[ schedule,JacTU ] = controls2Schedule( uk,schedule,'uScale',opt.uScale,...
    'partials',opt.partials);



targetObj = callArroba(target,{state,wellSol,netSol,schedule},'ComputePartials', opt.partials,'leftSeed',opt.leftSeed);
    

f = double(targetObj);
Jac = [];
if opt.partials
    [nSG] = nGridStateVariables( system.activeComponents );

	Jx = cell2mat(targetObj.jac(1:nSG))*JacTX;
    Jv = [cell2mat(targetObj.jac(nSG+1:2*nSG+1))*JacTW,cell2mat(targetObj.jac(2*nSG+2))*JacTN];
    Ju = cell2mat(targetObj.jac(2*nSG+3))*cell2mat(JacTU);
    
    if ~(size(opt.xRightSeeds,1)==0)
        Jac.J = Jx*opt.xRightSeeds + Ju*opt.uRightSeeds + Jv*opt.vRightSeeds;
    else
        Jac.Jx = Jx;
        Jac.Jv = Jv;
        Jac.Ju = Ju;
    end        
end


end