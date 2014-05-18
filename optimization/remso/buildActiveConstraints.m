function [ As,bs ] = buildActiveConstraints(A,b,active,sign)
%
%  As = sign*A(active,:)
%  bs = sign*b(active)
%
guardx =@(z)ifEmptyZero(z,size(A{1,1}));

ActiveMatrixInput = repmat(active,1,size(A,2));

As = cellfun(@(a1,a2)sign*subsref(guardx(a1),struct('type','()','subs',{{a2,':'}})),...
    A,ActiveMatrixInput,'UniformOutput',false);
bs = cellfun(@(a1,a2)sign*subsref(a1,struct('type','()','subs',{{a2}})),...
    b,active,'UniformOutput',false);


end

