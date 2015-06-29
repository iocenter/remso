function [ ok,callCount,k] = checkSendReceive(nRanks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ok = true;
k = 0;


for i = 0:nRanks-1
    for j = setdiff(0:nRanks-1,i)
        
        callCount = 1;
        
        while exist([num2str(i),'_',num2str(j),'_',num2str(callCount),'_Send_MPI.mat']) && ok
                k = k+1;
                loadmat = load([num2str(i),'_',num2str(j),'_',num2str(callCount),'_Send_MPI']);
                mat = loadmat.mat;
                varsSend = mat;
                
                loadmat = load([num2str(i),'_',num2str(j),'_',num2str(callCount),'_Recv_MPI']);
                mat = loadmat.mat;
                varsRecv = mat;
            
                    
            n = numel(varsSend);
            masterVal = reshape(varsSend,n,1) ;
            
                ok = ok & all(masterVal == reshape(varsRecv,n,1));
                if ~ok
                    break;
                end
            
            callCount = callCount +1 ;
        end
              
        callCount = callCount - 1 ;
        if ~ok
            break;
        end
    end
    if ~ok
        break;
    end
end

end

