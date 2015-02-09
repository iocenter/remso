
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
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/testFunctions'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('tankModel'));




predictionTime = 600;
totalPredictionSteps = 100;
totalControlSteps = 50;

simPerCtrl = totalPredictionSteps/totalControlSteps;

stepControl = cell2mat(arrayfun(@(index)ones(simPerCtrl,1)*index,(1:totalControlSteps)','UniformOutput',false));

ci = @(kk)controlIncidence(stepControl,kk);

dt = predictionTime/totalPredictionSteps;

state = [ 0 , 0.4, 1 ]';
u1 = 0.01;
u = repmat({u1},totalControlSteps,1);


ss.step = repmat({@(xS,u,varargin) tankAcadoModelAlg(xS,u,dt,varargin{:})},totalPredictionSteps,1);
ss.ci = ci;
ss.nv = 1;
ss.state = state;

obj = cell([],totalPredictionSteps,1);
for k =  1:totalPredictionSteps
    obj{k} = @(x,u,v,varargin) objectiveTest2(k,x,u,v,totalPredictionSteps,varargin{:});
end
targetObj = @(x,u,v,varargin) sepTarget(x,u,v,obj,ss,varargin{:});



maxError = unitTest(u,ss,obj,'totalSteps',10)



