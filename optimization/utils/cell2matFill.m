function [ matrix ] = cell2matFill( cellMatrix,iDims,jDims)
% similar to cell2mat but fill the empty cells with zeros if the given size
%
%

[di,dj] = size(cellMatrix);

if nargin < 2
    celldim = size(cellMatrix{1,1});
    iDims = celldim(1)*ones(di,1);
    jDims = celldim(2)*ones(1,dj);
end

matrix = zeros(sum(iDims),sum(jDims));


for i = 1:di
    for j = 1:dj
        if ~isempty(cellMatrix{i,j})
            matrix(sum([1;iDims(1:i-1)]):sum([0;iDims(1:i)]),sum([1,jDims(1:j-1)]):sum([0,jDims(1:j)])) = cellMatrix{i,j};
        end
    end
end
matrix = sparse(matrix);



end

