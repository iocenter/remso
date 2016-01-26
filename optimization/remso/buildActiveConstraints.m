function [ As,bs,rs ] = buildActiveConstraints(A,b,active,sign,varargin)
%
%  As = sign*A(active,:)
%  bs = sign*b(active)
%  rs = sign*R(active)
%
opt = struct('R',[]);
opt = merge_options(opt, varargin{:});

rs = [];

[m,n] = size(A);

dimZ = repmat(max(cellfun(@(x)size(x,1),A),[],2),1,n);
dimU = repmat(max(cellfun(@(x)size(x,2),A),[],1),m,1);

dim = arrayfun(@(dz,du)[dz,du],dimZ,dimU,'UniformOutput',false);


guardx =@(z,d)ifEmptyZero(z,d);

ActiveMatrixInput = repmat(active,1,size(A,2));

As = cellfun(@(a1,a2,d)sign*subsref(guardx(a1,d),struct('type','()','subs',{{a2,':'}})),...
    A,ActiveMatrixInput,dim,'UniformOutput',false);
bs = cellfun(@(a1,a2,d)sign*subsref(a1,struct('type','()','subs',{{a2}})),...
    b,active,'UniformOutput',false);

if ~isempty(opt.R)
    rs = cellfun(@(a1,a2)sign*subsref(a1,struct('type','()','subs',{{a2}})),...
        opt.R,active,'UniformOutput',false);
end


end

