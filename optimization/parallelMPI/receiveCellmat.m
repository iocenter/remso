function [ mat ] = receiveCellmat(nC,source)

%receivedimension
ijDims = NMPI_Recv(nC+1,source);
iDims = ijDims(1:nC);
jDim = ijDims(end);

iSum = sum(iDims);
n = iSum*jDim;
mat = NMPI_Recv(n,source);
mat = reshape(mat,iSum,jDim);
mat = mat2cell(mat,iDims,jDim);




end

