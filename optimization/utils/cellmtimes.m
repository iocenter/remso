function [ y ] = cellmtimes( A,b,varargin )
%
%  y = A * b
%
opt = struct('lowerTriangular',false,'ci',[]);
opt = merge_options(opt, varargin{:});

if isempty(opt.ci)
    opt.ci = @(kk)controlIncidence([],kk);
end

[mc1,nc1] = size(A);

mic = size(A{1,1},1);
bdim= size(b{1},2);
y = repmat({zeros(mic,bdim)},mc1,1);

if opt.lowerTriangular
    for ic = 1:mc1
        mjc = callArroba(opt.ci,{ic});
        for jc=1:mjc
            y{ic} = y{ic} + A{ic,jc}*b{jc};
        end
    end
else
    for ic = 1:mc1
        for jc=1:nc1
            y{ic} = y{ic} + A{ic,jc}*b{jc};
        end
    end
end