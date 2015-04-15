function [ obj,objJac] = objSumRisks(s,u,varargin)

opt = struct('gradients',false,'leftSeed',[],'sRightSeed',[],'uRightSeed',[]);
opt = merge_options(opt, varargin{:});


obj = sum(s);

objJac = [];
if opt.gradients
    uDims = cellfun(@numel,u);

    if size(opt.leftSeed,2) == 0 && size(opt.sRightSeed,1) == 0

        objJac.Js = ones(1,numel(s));
        objJac.Ju = mat2cell(sparse(1,sum(uDims)),1,uDims);
        
    elseif size(opt.leftSeed,2) ~= 0 && size(opt.sRightSeed,1) == 0
    
        objJac.Js = repmat(opt.leftSeed,1,numel(s));
                 %= opt.leftSeed*ones(1,numel(s))
        objJac.Ju = mat2cell(sparse(size(opt.leftSeed,1),sum(uDims)),1,uDims);

                 
    elseif size(opt.leftSeed,2) == 0 && size(opt.sRightSeed,1) ~= 0
        
        objJac.J = sum(opt.sRightSeed,2);
                %= ones(1,numel(s))*opt.sRightSeed + 0*opt.uRightSeed
    else
        error('not supported');
    end
end