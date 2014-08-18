function [f,Jac] = mrstTimePointFuncWrapper(xfk,uk,vk,target,schedule,wellSol,varargin)
%
% Interface with Remso a mrstPointFunction
%
% SYNOPSIS:
%  [f,Jac] = mrstTimePointFuncWrapper(xfk,uk,vk,target,schedule,wellSol,)
%  [f,Jac] = mrstTimePointFuncWrapper(xfk,uk,vk,target,schedule,wellSol, 'pn', pv, ...)
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
%   Jac - Jacobain of the targe function
%
% SEE ALSO:
%
%

opt = struct('partials',false,'leftSeed',[],'xScale',[],'vScale',[],'uScale',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});

[wellSol] = algVar2wellSol(vk,wellSol,'vScale',opt.vScale);
[ state,JacTX] = stateVector2stateMrst( xfk,'xScale',opt.xScale,...
    'partials',opt.partials);
[ schedule,JacTU ] = controls2Schedule( uk,schedule,'uScale',opt.uScale,...
    'partials',opt.partials);


targetObj = callArroba(target,{state,wellSol,schedule},'ComputePartials', opt.partials,'leftSeed',opt.leftSeed);
    

f = double(targetObj);
Jac = [];
if opt.partials
    Jv = bsxfun(@times,cell2mat(targetObj.jac(3:5)),opt.vScale');
    Jx = cell2mat(targetObj.jac(1:2))*JacTX;
    Ju = cell2mat(targetObj.jac(6))*JacTU;
    
    if ~isempty(opt.xRightSeeds)
        Jac.J = Jx*opt.xRightSeeds + Ju*opt.uRightSeeds + Jv*opt.vRightSeeds;
    else
        Jac.Jx = Jx;
        Jac.Jv = Jv;
        Jac.Ju = Ju;
    end        
end


end