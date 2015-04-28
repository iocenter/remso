function [o,Jac] = npvStages(x,u,v,varargin)
% we assume the (negative) of the cashflow (NPV) being accumulated in the last algebraic

opt = struct('gradients',false,'leftSeed',[],'vRightSeeds',[],'xRightSeeds',[],'uRightSeeds',[]);
opt = merge_options(opt, varargin{:});


vDims = cellfun(@numel,v);
xDims = cellfun(@numel,x);
uDims = cellfun(@numel,u);

o = cellfun(@(vi)vi(end),v);

Jac = [];
if opt.gradients
    Jv = cellfun(@(vi,k)sparse(k,numel(vi),1,numel(v),numel(vi)),v',num2cell(1:numel(v)), 'UniformOutput',false);  %% this is the jacobian of the output function
    if size(opt.leftSeed,2)~=0

        Jac.Jv = mat2cell(opt.leftSeed*cell2mat(Jv),size(opt.leftSeed,1),vDims);
        Jac.Jx = mat2cell(sparse(size(opt.leftSeed,1),sum(xDims)),size(opt.leftSeed,1),xDims);
        Jac.Ju = mat2cell(sparse(size(opt.leftSeed,1),sum(uDims)),size(opt.leftSeed,1),uDims);

    elseif iscell(opt.vRightSeed) && size(opt.vRightSeed{1},1)~=0

        Jac.J = cell2mat(Jv)*cell2mat(opt.vRightSeed);
        
    else
        Jac.Jv = Jv;
        Jac.Jx = mat2cell(sparse(sum(vDims),sum(xDims)),sum(vDims),xDims);
        Jac.Ju = mat2cell(sparse(sum(vDims),sum(uDims)),sum(vDims),uDims);
    end
end


end

