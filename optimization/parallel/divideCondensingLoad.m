function [workerCondensingScheduleW,clientCondensingSchedule,uStart,workerLoad,avgW] = divideCondensingLoad(totalPredictionSteps,ci,uDims,nWorkers)
% This can be done much better, try different strategies if needed

%  Returns a job-division for condensing per-worker
%  f
%
%

totalControlSteps = length(uDims);
workerCondensingSchedule = cell(1,nWorkers);
workerLoad = zeros(1,nWorkers);


[ uLoad,uStart ] = computeControlLoad(ci,totalControlSteps,totalPredictionSteps);
uGroupLoad = uLoad.*uDims';

avgW = sum(sum(uGroupLoad)+totalPredictionSteps)/nWorkers;


clientCondensingSchedule.correction = 1;
workerCondensingSchedule{1} = [workerCondensingSchedule{1},0]; % correction calculation
workerLoad(1) = totalPredictionSteps;

clientCondensingSchedule.control = zeros(1,totalControlSteps);

for taskGroup = 1:totalControlSteps
   
    
    taskGroupLoad = uGroupLoad(taskGroup);
    
    
	[~,w] = min(workerLoad);

	workerLoad(w) = workerLoad(w) + taskGroupLoad;  

    
    workerCondensingSchedule{w} = [workerCondensingSchedule{w},taskGroup];
    clientCondensingSchedule.control(taskGroup) = w;
    
end


workerCondensingScheduleW = Composite();
for w = 1:nWorkers
    workerCondensingScheduleW{w} = workerCondensingSchedule{w};
end



end