function [ obj,objJac] = objSumRisks(s,varargin)

opt = struct('gradients',false,'leftSeed',[],'sRightSeed',[]);
opt = merge_options(opt, varargin{:});

obj = sum(s);

objJac = [];
if opt.gradients

    if size(opt.leftSeed,2) == 0 && size(opt.sRightSeed,1) == 0

        objJac.Js = ones(1,numel(s));
        
    elseif size(opt.leftSeed,2) ~= 0 && size(opt.sRightSeed,1) == 0
    
        objJac.Js = repmat(opt.leftSeed,1,numel(s));
                 %= opt.leftSeed*ones(1,numel(s))
                 
    elseif size(opt.leftSeed,2) == 0 && size(opt.sRightSeed,1) ~= 0
        
        objJac.J = sum(opt.sRightSeed,2);
                %= ones(1,numel(s))*opt.sRightSeed
    else
        error('not supported');
    end
end