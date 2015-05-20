function [ varDistributed ] = createEmptyDistributedVar( jobSchedule )

varDistributed = cell(numel(jobSchedule.work2Job{jobSchedule.my_rank+1}),1);

end






