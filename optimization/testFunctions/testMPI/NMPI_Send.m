function [ ] = NMPI_Send( mat,n,dest )

global callSendCount
callSendCount(dest+1) = callSendCount(dest+1) + 1;

assert(numel(mat) == n);

myRank = NMPI_Comm_rank();

name = [num2str(myRank),'_',num2str(dest),'_',num2str(callSendCount(dest+1)),'_Send_MPI'];
save(name,'mat');


path = fileparts(which('NMPI_Send'));
rmpath(path);
NMPI_Send(mat,n,dest);
addpath(path);



end

