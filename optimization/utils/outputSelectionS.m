function [ obj,objJac] = outputSelectionS(s,u,index,varargin)

opt = struct('gradients',false,'leftSeed',[],'sRightSeed',[],'uRightSeed',[]);
opt = merge_options(opt, varargin{:});


obj = s(index);

objJac = [];
if opt.gradients
    ns = sum(index);
    Js = sparse(1:sum(index),find(index),1,ns,numel(s));
   
    uDims = cellfun(@numel,u);

    if size(opt.leftSeed,2) == 0 && size(opt.sRightSeed,1) == 0

        objJac.Js = Js;
        objJac.Ju = mat2cell(sparse(ns,sum(uDims)),ns,uDims);
        
    elseif size(opt.leftSeed,2) ~= 0 && size(opt.sRightSeed,1) == 0
    
        objJac.Js = opt.leftSeed*Js;
        objJac.Ju = mat2cell(sparse(size(opt.leftSeed,1),sum(uDims)),size(opt.leftSeed,1),uDims);

                 
    elseif size(opt.leftSeed,2) == 0 && size(opt.sRightSeed,1) ~= 0
        
        objJac.J = Js*opt.sRightSeed;
    else
        error('not supported');
    end
end
