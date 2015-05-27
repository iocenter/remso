function [ ok,callCount,k] = checkSendReceive(nRanks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ok = true;

varsSend = cell(nRanks,1);
varsRecv = cell(nRanks,1);


for i = 1:nRanks
    for j = setdiff(1:nRanks,i)
        
        callCount = 1;
        
        while exist([num2str(i),'_',num2str(j),'_',num2str(callCount),'_Send_MPI.mat']) && ok
            for k = 0:nRanks-1
                loadmat = load([num2str(i),'_',num2str(j),'_',num2str(callCount),'_Send_MPI']);
                mat = loadmat.mat;
                varsSend{k+1} = mat;
                
                loadmat = load([num2str(i),'_',num2str(j),'_',num2str(callCount),'_Recv_MPI']);
                mat = loadmat.mat;
                varsRecv{k+1} = mat;
            end
                    
            n = numel(varsSend{1});
            masterVal = reshape(varsSend{1},n,1) ;
            
            for k = 1:nRanks
                ok = ok & all(masterVal == reshape(varsRecv{k},n,1));
                if ~ok
                    break;
                end
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

