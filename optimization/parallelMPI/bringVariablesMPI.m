function [ var ] = bringVariablesMPI( distVar,jobSchedule  )
%
%  bring back variables to the client
%
%
num_ranks = jobSchedule.num_ranks;
my_rank = jobSchedule.my_rank;
Master_rank = jobSchedule.Master_rank;

var = cell(size(jobSchedule.job2Work,1),1);

if iscell(distVar)
    isCellCell = iscell(distVar{1});  %more syncro needed
    if my_rank == Master_rank
        for w = 1:num_ranks
            if w == Master_rank+1
                var(jobSchedule.work2Job{w}) = distVar;
            else
                nJobs = numel(jobSchedule.work2Job{w});
                
                mat = receiveCell(nJobs,w-1,isCellCell);
                
                var(jobSchedule.work2Job{w}) = mat;
                
            end
        end
    else
        
        sendCell(distVar,Master_rank,isCellCell);
        
    end
    
    
else
    error('notImplemented')
end




end

function [mat] = receiveCell(nJobs,rank,isCellCell)
if isCellCell
    iDims = NMPI_Recv(nJobs,rank);
    mat = receiveCellmat(iDims,rank);
    mat = mat2cell(mat,iDims,1);
else
    mat = receiveCellmat(nJobs,rank);
end

end

function [] = sendCell(distVar,Master_rank,isCellCell)
if isCellCell
    cDims = cellfun(@numel,distVar);
    NMPI_Send(cDims,numel(cDims),Master_rank);
    distVar = vertcat(distVar{:});
end
sendCellmat(distVar,Master_rank);


end

