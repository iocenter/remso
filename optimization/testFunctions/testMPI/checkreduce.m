function [ ok,callCount,k] = checkreduce(nRanks )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ok = true;

callCount = 1;

varsSend = cell(nRanks,1);
varsRecv = cell(nRanks,1);

opSend = cell(nRanks,1);
opRecv = cell(nRanks,1);

while exist([num2str(0),'_',num2str(callCount),'_red_sent','MPI.mat']) && ok
    for k = 0:nRanks-1
        loadmat = load([num2str(k),'_',num2str(callCount),'_red_sent','MPI']);
        mat = loadmat.mat;
        varsSend{k+1} = mat;
        opSend{k+1} = loadmat.op;
        
        loadmat = load([num2str(k),'_',num2str(callCount),'_red_recv','MPI']);
        mat = loadmat.mat;
        varsRecv{k+1} = mat;
        opRecv{k+1} = loadmat.op;
    end
    
    n = numel(varsSend{1});
    masterVal = reshape(varsSend{1},n,1) ;
    op = opSend{1};
    
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
    
    val = masterVal;
    z = zeros(n,1);
    for k = 2:nRanks
        val = opF(val,varsSend{k});
        ok = ok && all(z == reshape(varsRecv{k},n,1));
    end
    ok = ok && all(val == reshape(varsRecv{1},n,1));
    
    for k = 1:nRanks;
        ok = ok && strcmp(op,opSend{k}) && strcmp(op,opRecv{k});
    end
    
    
    callCount = callCount +1 ;
end
callCount = callCount - 1 ;


end

