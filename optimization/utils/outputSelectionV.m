function [o,Jac] = outputSelectionV(x,u,v,selection,varargin)

opt = struct('gradients',false,'leftSeed',[],'vRightSeeds',[],'xRightSeeds',[],'uRightSeeds',[]);
opt = merge_options(opt, varargin{:});


o = cell2mat(cellfun(@(vi)vi(selection),v,'UniformOutput',false));

Jac = [];
if opt.gradients
    
    vDims = cellfun(@numel,v);
    xDims = cellfun(@numel,x);
    uDims = cellfun(@numel,u);

    nS = sum(selection);
    nv = numel(v);
    no = nv*nS;
    
    Jv = arrayfun(@(vDim,i)sparse((i-1)*nS+1:i*nS,find(selection),1,no,vDim),vDims',1:nv, 'UniformOutput',false);  %% this is the jacobian of the output function
    if size(opt.leftSeed,2)~=0

        Jac.Jv = mat2cell(opt.leftSeed*cell2mat(Jv),size(opt.leftSeed,1),vDims);
        Jac.Jx = mat2cell(sparse(size(opt.leftSeed,1),sum(xDims)),size(opt.leftSeed,1),xDims);
        Jac.Ju = mat2cell(sparse(size(opt.leftSeed,1),sum(uDims)),size(opt.leftSeed,1),uDims);

    elseif iscell(opt.vRightSeeds) && size(opt.vRightSeeds{1},1)~=0
        
        Jac.J = cell2mat(Jv)*cell2mat(opt.vRightSeeds);
        
    else
        Jac.Jv = Jv;
        Jac.Jx = mat2cell(sparse(no,sum(xDims)),no,xDims);
        Jac.Ju = mat2cell(sparse(no,sum(uDims)),no,uDims);
    end
end


end

