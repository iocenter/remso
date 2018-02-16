function [f,Jac] = zeroTarget(x,u,v,varargin)

opt = struct('gradients',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[],'usliced',[]);
opt = merge_options(opt, varargin{:});

if isempty(x)
    x = {};
end
if isempty(v)
    v = {};
end
if isempty(u)
    u = {};
end

xDims = cellfun(@numel,x);
vDims = cellfun(@numel,v);
uDims = cellfun(@numel,u);

f =  0;

Jac = [];
if opt.gradients
    if ~isempty(opt.xRightSeeds) && ~isempty(opt.xRightSeeds{1})
        Jac.J = zeros(size(f,1),size(opt.xRightSeeds{1},2));
    else
        nL = size(opt.leftSeed,2);
        
        Jac.Ju = mat2cell(zeros(nL,sum(uDims)),nL,uDims );
        Jac.Jx = mat2cell(zeros(nL,sum(xDims)),nL,xDims );
        Jac.Jv = mat2cell(zeros(nL,sum(vDims)),nL,vDims );
       
    end
end

end