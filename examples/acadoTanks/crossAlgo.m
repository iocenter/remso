
clc
clear
clear global

% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/parallel'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoCrossSequential'));
addpath(genpath('../../optimization/remsoCross'));
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('tankModel'));


predictionTime = 600;
totalPredictionSteps = 100;
totalControlSteps = 50;

simPerCtrl = totalPredictionSteps/totalControlSteps;

stepControl = cell2mat(arrayfun(@(index)ones(simPerCtrl,1)*index,(1:totalControlSteps)','UniformOutput',false));

ci = arroba(@controlIncidence,2,{stepControl});

dt = predictionTime/totalPredictionSteps;



initPool()

%%  Who will do what - Distribute the computational effort!
nWorkers = getNumWorkers();
if nWorkers == 0
    nWorkers = 1;
end
[ jobSchedule ] = divideJobsSequentially(totalPredictionSteps ,nWorkers);
jobSchedule.nW = nWorkers;

work2Job = Composite();
for w = 1:nWorkers
    work2Job{w} = jobSchedule.work2Job{w};
end



stepClient = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    ss.stepClient{k} = @(x0,u,varargin) tankAcadoModelAlg(x0,u,dt,varargin{:});
end


spmd
    
    stepW = cell(totalPredictionSteps,1);
    for k=1:totalPredictionSteps
        cik = callArroba(ci,{k});
        stepW{k} = arroba(@tankAcadoModelAlg,...
            [1,2],...
            {dt},...
            true);
    end
    step = stepW;
end

ss.state = [ 0 , 0.4, 1 ]';
ss.nv = 1;
ss.jobSchedule = jobSchedule;
ss.work2Job = work2Job;
ss.step = step;
ss.ci = ci;


u1 = 0.01;




% the objective function is a separable, exploit this!
spmd
    nJobsW = numel(work2Job);
    objW = cell(nJobsW,1);
    for i = 1:nJobsW
        k = work2Job(i);
        cik = callArroba(ci,{k});
        
        objW{i} = arroba(@objectiveTest2,...
            [2,3,4],...
            {...
            k,...
            totalPredictionSteps...
            },true);
    end
    obj = objW;
end
targetObj = @(x,u,v,varargin) sepTarget(x,u,v,obj,ss,jobSchedule,work2Job,varargin{:});


%%% objective function on the client side (for plotting!)
objClient = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    objClient{k} = arroba(@objectiveTest2,...
        [2,3,4],...
        {...
        k,...
        totalPredictionSteps...
        },true);
end



lbx = repmat({[0;0.1;0.1]},totalPredictionSteps,1);
ubx = repmat({[inf;5;5]},totalPredictionSteps,1);

lbu = repmat({0},totalControlSteps,1);
ubu = repmat({0.3},totalControlSteps,1);

u = repmat({u1},totalControlSteps,1);
plotFunc = @(xi,ui,v,di,varargin) plotSolution(xi,ui,ss,objClient,'simulate',false,varargin{:});


[u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbu',lbu,'ubu',ubu,'lkMax',10,'plotFunc',plotFunc,'max_iter',200,'tol',1e-7,'plot',false,'lkMax',20);


plotFunc(x,u,[],xd,'simulate',true)

