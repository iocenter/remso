function [obj,Jac] = objectiveTest2( k,x,u,v,totalSteps,varargin)
%objecitive being accumulated in the first state
%  fV(k),dfdxk,dfdui,  dfdvk(x{k},ui,v{k},'partials',true,

opt = struct('partials',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});

notEmptyV = ~isempty(v);
if notEmptyV
    obj =(x(2)-2)^2 + u^2 + v^2*0;
else
    obj =(x(2)-2)^2 + u^2; 
end
Jx = [0,2*(x(2)-2),0];
Ju = 2*u;
if notEmptyV
    Jv = 2*v*0;
end

if ~(size(opt.xRightSeeds,1)==0)
    if notEmptyV
        Jac.J = Jx*opt.xRightSeeds + Ju*opt.uRightSeeds + Jv*opt.vRightSeeds ;
    else
        Jac.J = Jx*opt.xRightSeeds + Ju*opt.uRightSeeds  ;
    end
    
elseif ~(size(opt.leftSeed,2)==0)
    Jac.Jx = opt.leftSeed*Jx;
    Jac.Ju = opt.leftSeed*Ju;
    if notEmptyV
        Jac.Jv = opt.leftSeed*Jv;
    else
        Jac.Jv = [];
    end
else
    Jac.Jx = Jx;
    Jac.Ju = Ju;
    if notEmptyV
        Jac.Jv = Jv;
    else
        Jac.Jv = [];
    end
end



end

