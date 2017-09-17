
clc
clear
clear global

here = fileparts(mfilename('fullpath'));
if isempty(here)
    here = pwd();
end



% Include REMSO functionalities
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'multipleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remso')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoCrossSequential')));
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


ss.step = repmat({@(xS,u,varargin) tankAcadoModelAlg(xS,u,dt,varargin{:})},totalPredictionSteps,1);
ss.ci = ci;
ss.state = state;

obj = cell([],totalPredictionSteps,1);
for k =  1:totalPredictionSteps
    obj{k} = @(x,u,v,varargin) objectiveTest2(k,x,u,v,totalPredictionSteps,'scale',0.1,varargin{:});
end
targetObj = @(x,u,v,varargin) sepTarget(x,u,v,obj,ss,varargin{:});


lbx = repmat({[0;0.1;0.1]},totalPredictionSteps,1);
ubx = repmat({[inf;5;5]},totalPredictionSteps,1);

lbu = repmat({0.01},totalControlSteps,1);
ubu = repmat({0.3},totalControlSteps,1);

u = repmat({u1},totalControlSteps,1);
plotFunc = @(xi,ui,v,di,varargin) plotSolution(xi,ui,ss,obj,'simulate',false,varargin{:});


[u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbu',lbu,'ubu',ubu,'lkMax',10,'plotFunc',plotFunc,'max_iter',200,'tol',1e-5,'plot',false,'lkMax',20,...
                              'saveIt',true,...
                              'lbxH',lbx,'ubxH',ubx,'condense',false);


plotFunc(x,u,[],xd,'simulate',true)

