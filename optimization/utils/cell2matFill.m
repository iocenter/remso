function [ matrix ] = cell2matFill( cellMatrix,celldim)
% similar to cell2mat but fill the empty cells with zeros if the given size
%
%

if nargin < 2
    celldim = size(cellMatrix{1,1});
end

[di,dj] = size(cellMatrix);

matrix = zeros(di*celldim(1),dj*celldim(2));


for i = 1:di
    for j = 1:dj
        if ~isempty(cellMatrix{i,j})
            matrix((i-1)*celldim(1)+1:(i)*celldim(1),(j-1)*celldim(2)+1:(j)*celldim(2)) = cellMatrix{i,j};
        end
    end
end
matrix = sparse(matrix);



end

