function [ mat ] = NMPI_Recv( lengthMat,source )

global callRecvCount
callRecvCount(source+1) = callRecvCount(source+1) + 1;

myRank = NMPI_Comm_rank();

name = [num2str(source),'_',num2str(myRank),'_',num2str(callRecvCount(source+1)),'_Recv_MPI'];

path = fileparts(which('NMPI_Recv'));
rmpath(path);
mat = NMPI_Recv(lengthMat,source);
addpath(path);

save(name,'mat');

assert(numel(mat) == lengthMat);


end

