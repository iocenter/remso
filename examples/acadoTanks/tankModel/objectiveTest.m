function [obj,Jac] = objectiveTest( k,x,u,v,totalSteps,varargin)
%objecitive being accumulated in the first state

opt = struct('partials',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});

notEmptyV = ~isempty(v);

if(totalSteps == k)
    obj = x(1);%*0 + (u-0.02)^2;
    Jx = zeros(1,numel(x));
    Ju = zeros(1,numel(u));%+2*(u-0.02);
    Jx(1) = 1;%*0;
else
    obj = 0;%+(u-0.001)^2;
    Jx = zeros(1,numel(x));
    Ju = zeros(1,numel(u));%+2*(u-0.02);
end


if notEmptyV
    Jv = 2*v*0;
end

if ~isempty(opt.xRightSeeds)
    if notEmptyV
        Jac.J = Jx*opt.xRightSeeds + Ju*opt.uRightSeeds + Jv*opt.vRightSeeds ;
    else
        Jac.J = Jx*opt.xRightSeeds + Ju*opt.uRightSeeds  ;
    end
    
elseif ~isempty(opt.leftSeed)
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

