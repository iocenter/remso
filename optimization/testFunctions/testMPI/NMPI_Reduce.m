function [ mat ] = NMPI_Reduce(mat,lengthMat,op,masterRank )



global callReduceCount
callReduceCount = callReduceCount + 1;

myRank = NMPI_Comm_rank();

assert(numel(mat) == lengthMat)

name = [num2str(myRank),'_',num2str(callReduceCount),'_red_sent','MPI'];
save(name,'mat','op');

path = fileparts(which('NMPI_Reduce'));
rmpath(path);
mat = NMPI_Reduce(mat,lengthMat,op,masterRank);
addpath(path);

name = [num2str(myRank),'_',num2str(callReduceCount),'_red_recv','MPI'];
save(name,'mat','op');




end

