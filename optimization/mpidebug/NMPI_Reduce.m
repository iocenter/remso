function [ opAllMat ] = NMPI_Reduce(mat,lengthMat,op,masterRank )

nRanks = NMPI_Comm_size();
myRank = NMPI_Comm_rank();

if masterRank == myRank 
    
    allMats = cell(nRanks,1);
    for r = 0:nRanks-1
        if r ~= myRank
            allMats{r+1} = NMPI_Recv(lengthMat,r);
        else
            allMats{r+1} = mat;
        end
    end
      
    switch op
        case {'+'}
            opF = @plus;
        case {'*'}
            opF = @times;
        case {'N'}
            opF = @min;
        case {'M'}
            opF = @max;
        otherwise
            error('Unknown op.')
    end
    
    opAllMat = allMats{1};
    for k = 2:numel(allMats)
        opAllMat = opF(opAllMat,allMats{k});
    end
       
else 
    NMPI_Send(mat,lengthMat,masterRank);
    opAllMat = zeros(size(mat)); 
end



end

