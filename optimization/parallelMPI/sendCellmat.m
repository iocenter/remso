function [] = sendCellmat( cellmat,dest )


iDims = cellfun(@(x)size(x,1),cellmat);
jDim = size(cellmat{1},2);

lengthI = numel(iDims);
NMPI_Send([iDims;jDim],sum(lengthI)+1,dest);

mat = cell2mat(cellmat);
n = numel(mat);
mat = reshape(mat,n,1);
NMPI_Send(mat,n,dest);


end

