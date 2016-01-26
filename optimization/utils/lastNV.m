function [o,Jac] = lastNV(x,u,v,n,varargin)

opt = struct('gradients',false,'leftSeed',[],'vRightSeeds',[],'xRightSeeds',[],'uRightSeeds',[]);
opt = merge_options(opt, varargin{:});


o = cell2mat(cellfun(@(vi)vi(end-n+1:end),v,'UniformOutput',false));

Jac = [];
if opt.gradients
    
    vDims = cellfun(@numel,v);
    xDims = cellfun(@numel,x);
    uDims = cellfun(@numel,u);

    nv = numel(v);
    no = nv*n;
    
    Jv = arrayfun(@(vDim,i)sparse((i-1)*n+1:i*n,vDim-n+1:vDim,1,no,vDim),vDims',1:nv, 'UniformOutput',false);  %% this is the jacobian of the output function
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

