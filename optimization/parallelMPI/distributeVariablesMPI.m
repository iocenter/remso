function [ varDistributed ] = distributeVariablesMPI( var,jobSchedule )
%DISTRIBUTEVARIABLESMPI Summary of this function goes here
%   Detailed explanation goes here


my_rank = jobSchedule.my_rank;
Master_rank = jobSchedule.Master_rank;


if my_rank == Master_rank
    for w = 1:numel(jobSchedule.work2Job)
        if my_rank + 1  == w
            varDistributed = var(jobSchedule.work2Job{w});
        else
            sendCellmat(var(jobSchedule.work2Job{w}),w-1);
        end
    end
else
    nJobs = numel(jobSchedule.work2Job{my_rank+1});
    varDistributed = receiveCellmat(nJobs,Master_rank);
end


