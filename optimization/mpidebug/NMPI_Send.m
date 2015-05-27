function [ ] = NMPI_Send( mat,n,dest )

assert(numel(mat) == n);

myRank = NMPI_Comm_rank();

name = [num2str(myRank),'_',num2str(dest),'_','MPI'];
nameAck = [name '.ack'];

save(name,'mat');

fid = fopen(nameAck,'w');
fclose(fid);

block = fopen(nameAck,'r');
while block ~= -1
    fclose(block);
    pause(0.01);
    block = fopen(nameAck,'r');
end

delete([name '.mat']);


end

