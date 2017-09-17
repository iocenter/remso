
clc
clear
clear global


here = fileparts(mfilename('fullpath'));
if isempty(here)
    here = pwd();
end


% Include REMSO functionalities
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstDerivated')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'procedural')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'multipleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'parallel')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remso')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'testFunctions')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'utils')));
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



[maxError,crossError] = unitTest(u,ss,obj,'totalSteps',10)



