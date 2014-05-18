function [x,xf,v,u] = getVarsfromSimulation(simVars,controlSchedules,varargin)
warning('forwardStates does not contain the first simulation step anymore!')

opt = struct('xScale',[],'vScale',[],'uScale',[]);
opt = merge_options(opt, varargin{:});


totalPredictionSteps = numel(simVars);
totalControlSteps = numel(controlSchedules);


x = cell(totalPredictionSteps,1);
xf = cell(totalPredictionSteps,1);
v = cell(totalPredictionSteps,1);
u = cell(totalControlSteps,1);


for k =1:totalPredictionSteps
    xf{k} = stateMrst2stateVector(simVars{k}.forwardStates{2},'doScale',true,'xScale',opt.xScale );
end
for k =1:totalPredictionSteps-1
    x{k} = stateMrst2stateVector(simVars{k+1}.forwardStates{1},'doScale',true,'xScale',opt.xScale );
end
% by convention: The initial step in the last period+1 = last step
x{totalPredictionSteps} = stateMrst2stateVector(simVars{k+1}.forwardStates{1},'doScale',true,'xScale',opt.xScale );
for k =1:totalPredictionSteps
    v{k} = wellSol2algVar(simVars{k}.forwardStates{2}.wellSol,'doScale',true,'vScale',opt.vScale);
end


firstStep = 1;
for k = 1:totalControlSteps
    
    sch = wellSol2Schedule(simVars{firstStep}.forwardStates{2}.wellSol,controlSchedules(k));
    
    u{k} = schedule2Controls(sch,'doScale',true,'uScale',opt.uScale{k}); 
    
    firstStep = firstStep + numel(controlSchedules(k).step.val);
    
end


end

