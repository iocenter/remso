function [ y] = toStructuredCells( yVector,ny,varargin)
% transform a vector to a cell arrar divided containing ny elements each
% cell

opt = struct('T',false);
opt = merge_options(opt, varargin{:});

if ~opt.T
    y = mat2cell(yVector,ones(1,size(yVector,1)/ny)*ny,ones(1,size(yVector,2)));
else
    y = mat2cell(yVector,ones(1,size(yVector,1)),ones(1,size(yVector,2)/ny)*ny);
end

end

