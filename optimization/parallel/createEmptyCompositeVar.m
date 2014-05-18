function [ varDistributed ] = createEmptyCompositeVar(jobSchedule)
%
%  Create an empty distributed variable according to JobSchedule dimensions
%
%
varDistributed = Composite();
for wi = 1:length( varDistributed )
    nj = length(jobSchedule.work2Job{wi});
    varDistributed{wi} = cell(nj,1);
end

end

