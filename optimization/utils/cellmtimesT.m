function [ y ] = cellmtimesT( x,A,varargin )
%  y = x'A, exploit the structure of the predictor matrix
%  columnVector == true then x is a column vector, otherwise it is a line

opt = struct('lowerTriangular',false,'ci',[],'columnVector',true);
opt = merge_options(opt, varargin{:});

if isempty(opt.ci)
    opt.ci = @(kk)controlIncidence([],kk);
end

[mA,nA] = size(A);

nC = size(A{1,1},2);
y = repmat({zeros(1,nC)},1,nA);

if opt.columnVector
    if opt.lowerTriangular
        for iA = 1:mA
            maxjA = callArroba(opt.ci,{iA});
            for jA = 1:maxjA   %jA = jy
                y{jA} = y{jA} + x{iA}'*A{iA,jA};
            end
        end
    else
        for jA = 1:nA   %jA = jy
            for iA = 1:mA
                y{jA} = y{jA} + x{iA}'*A{iA,jA};
            end
        end
    end
else
    if opt.lowerTriangular
        for iA = 1:mA
            maxjA = callArroba(opt.ci,{iA});
            for jA = 1:maxjA   %jA = jy
                y{jA} = y{jA} + x{iA}*A{iA,jA};
            end
        end
    else
        for jA = 1:nA   %jA = jy
            for iA = 1:mA
                y{jA} = y{jA} + x{iA}*A{iA,jA};
            end
        end
    end
    
    
end

