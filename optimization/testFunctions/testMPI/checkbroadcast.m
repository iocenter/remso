function [ ok,callCount,k] = checkbroadcast(nRanks )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ok = true;

callCount = 1;

varsSend = cell(nRanks,1);
varsRecv = cell(nRanks,1);

while exist([num2str(0),'_',num2str(callCount),'_Bcast_sent','MPI.mat']) && ok
for k = 0:nRanks-1
    loadmat = load([num2str(k),'_',num2str(callCount),'_Bcast_sent','MPI']);
    mat = loadmat.mat;
    varsSend{k+1} = mat;
    
    loadmat = load([num2str(k),'_',num2str(callCount),'_Bcast_recv','MPI']);
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


end

