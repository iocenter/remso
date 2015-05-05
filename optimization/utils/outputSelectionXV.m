function [o,Jac] = outputSelectionXV(x,u,v,xS,vS,varargin)

opt = struct('gradients',false,'leftSeed',[],'vRightSeeds',[],'xRightSeeds',[],'uRightSeeds',[]);
opt = merge_options(opt, varargin{:});

o = cell2mat([cellfun(@(zi)zi(xS),x,'UniformOutput',false);
              cellfun(@(zi)zi(vS),v,'UniformOutput',false)]);

Jac = [];
if opt.gradients
    
    vDims = cellfun(@numel,v);
    xDims = cellfun(@numel,x);
    uDims = cellfun(@numel,u);
    
    nxS = sum(xS);
    nvS = sum(vS);

    nx = numel(x);
    nv = numel(v);
    nox = nx*nxS;
    no =  nox+nv*nvS;
    
    Jx = arrayfun(@(zDim,i)sparse(    (i-1)*nxS+1:    i*nxS,find(xS),1,no,zDim),xDims',1:nx, 'UniformOutput',false);
	Jv = arrayfun(@(zDim,i)sparse(nox+(i-1)*nvS+1:nox+i*nvS,find(vS),1,no,zDim),vDims',1:nv, 'UniformOutput',false);   

    if size(opt.leftSeed,2)~=0

        Jac.Jv = mat2cell(opt.leftSeed*cell2mat(Jv),size(opt.leftSeed,1),vDims);
        Jac.Jx = mat2cell(opt.leftSeed*cell2mat(Jx),size(opt.leftSeed,1),xDims);
        Jac.Ju = mat2cell(sparse(size(opt.leftSeed,1),sum(uDims)),size(opt.leftSeed,1),uDims);

    elseif iscell(opt.vRightSeeds) && size(opt.vRightSeeds{1},1)~=0
        
        Jac.J = cell2mat(Jx)*cell2mat(opt.xRightSeeds)+cell2mat(Jv)*cell2mat(opt.vRightSeeds);
        
    else
        Jac.Jv = Jv;
        Jac.Jx = Jx;
        Jac.Ju = mat2cell(sparse(no,sum(uDims)),no,uDims);
    end
end


end

