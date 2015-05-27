function [ mat] = NMPI_Bcast(mat,lengthMat,masterRank,myRank)

global callCount

callCount = callCount + 1;

assert(numel(mat) == lengthMat)
assert(myRank == NMPI_Comm_rank())


name = [num2str(myRank),'_',num2str(callCount),'_Bcast_sent','MPI'];
save(name,'mat');

path = fileparts(which('NMPI_Bcast'));
rmpath(path);
mat = NMPI_Bcast(mat,lengthMat,masterRank,myRank);
addpath(path);

name = [num2str(myRank),'_',num2str(callCount),'_Bcast_recv','MPI'];
save(name,'mat');



end
