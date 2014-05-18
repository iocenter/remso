function totalPredictionSteps = getTotalPredictionSteps(ss)

if isfield(ss,'jobSchedule')  %% we are executing the parallel version!
    totalPredictionSteps = size(ss.jobSchedule.job2Work,1);
else
    totalPredictionSteps = numel(ss.step);
end


end