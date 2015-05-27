function [ mat ] = NMPI_Recv( lengthMat,source )


myRank = NMPI_Comm_rank();

name = [num2str(source),'_',num2str(myRank),'_','MPI'];
nameAck = [name '.ack'];

block = fopen(nameAck,'r');
while block == -1
    pause(0.01);
    block = fopen(nameAck,'r');
end
fclose(block);

loadVars = load(name);
mat = loadVars.mat;

delete(nameAck)


assert(numel(mat) == lengthMat);


end

