function [ jobSchedule ] = divideJobsSequentially(totalPredictionSteps ,nWorkers )
%
%  Create a structure that encodes the data and job division of prediction
%  step computations in the number of available workers.
%
%  work2Job{w} gives a vector containing the indexes handled by worker w
%  job2Work(k,1) gives the worker in charge of k
%  job2Work(k,2) gives the index of work k in the corresponding worker


res = mod(totalPredictionSteps,nWorkers);
jobDiv = floor(totalPredictionSteps/nWorkers);

work2Job = cell(1,nWorkers);
job2Work = zeros(totalPredictionSteps,2);
for w = 1:nWorkers
   if w <= res
       work2Job{w} = (jobDiv+1)*(w-1)+1:(jobDiv+1)*w;
       
       job2Work((jobDiv+1)*(w-1)+1:(jobDiv+1)*w,1) = ones(jobDiv+1,1)*w;
       job2Work((jobDiv+1)*(w-1)+1:(jobDiv+1)*w,2) = 1:jobDiv+1;

   else
       work2Job{w} = (jobDiv)*(w-1)+1+res:(jobDiv)*w+res;
       
       job2Work((jobDiv)*(w-1)+1+res:(jobDiv)*w+res,1) =  ones(jobDiv,1)*w;
       job2Work((jobDiv)*(w-1)+1+res:(jobDiv)*w+res,2) = 1:jobDiv;
   end
   
end

jobSchedule.work2Job = work2Job;
jobSchedule.job2Work = job2Work;


end

