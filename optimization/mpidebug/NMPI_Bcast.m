function [ mat] = NMPI_Bcast(mat,lengthMat,masterRank,myRank)

assert(numel(mat) == lengthMat)

nRanks = NMPI_Comm_size();
myRank = NMPI_Comm_rank();


if masterRank == myRank % i'm sending
    for r = 0:nRanks-1
        if r ~= myRank
            NMPI_Send(mat,numel(mat),r)
        end
    end
else % i'm receiving
    mat = NMPI_Recv(lengthMat,masterRank);
end


end

