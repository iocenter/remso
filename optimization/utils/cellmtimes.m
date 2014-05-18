function [ y ] = cellmtimes( A,b,varargin )
%
%  y = A * b
%
opt = struct('lowerTriangular',false,'ci',@(kk)controlIncidence([],kk));
opt = merge_options(opt, varargin{:});

[mc1,nc1] = size(A);

mic = size(A{1,1},1);
y = repmat({zeros(mic,1)},mc1,1);

if opt.lowerTriangular
    for ic = 1:mc1
        mjc = opt.ci(ic);
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