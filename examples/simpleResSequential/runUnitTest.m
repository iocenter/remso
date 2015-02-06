% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global

% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi

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
addpath(genpath('reservoirData'));

%% Initialize reservoir -  the Simple reservoir
[reservoirP] = initReservoir( 'reallySimpleRes.data','Verbose',true);
%[reservoirP] = initReservoir( 'reallySimpleRes.data','Verbose',true);

% do not display reservoir simulation information!
mrstVerbose off;

% Number of reservoir grid-blocks
nCells = reservoirP.G.cells.num;

%% Multiple shooting problem set up
totalPredictionSteps = numel(reservoirP.schedule.step.val);  % MS intervals

% Schedule partition for each control period and for each simulated step
lastControlSteps = findControlFinalSteps( reservoirP.schedule.step.control );
controlSchedules = multipleSchedules(reservoirP.schedule,lastControlSteps);

stepSchedules = multipleSchedules(reservoirP.schedule,1:totalPredictionSteps);


% Piecewise linear control -- mapping the step index to the corresponding
% control 
ci  = arroba(@controlIncidence,2,{reservoirP.schedule.step.control});


%% Variables Scaling
xScale = setStateValues(struct('pressure',5*barsa,'sW',0.01),'nCells',nCells);


if (isfield(reservoirP.schedule.control,'W'))
    W =  reservoirP.schedule.control.W;
else
    W = processWellsLocal(reservoirP.G, reservoirP.rock,reservoirP.schedule.control(1),'DepthReorder', true);
end
wellSol = initWellSolLocal(W, reservoirP.state);
vScale = wellSol2algVar( wellSolScaling(wellSol,'bhp',5*barsa,'qWs',10*meter^3/day,'qOs',10*meter^3/day) );

cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
    'RATE',10*meter^3/day,...
    'ORAT',10*meter^3/day,...
    'WRAT',10*meter^3/day,...
    'LRAT',10*meter^3/day,...
    'RESV',0,...
    'BHP',5*barsa));


%% instantiate the objective function as an aditional Algebraic variable


%%% Last values is the objective
nCells = reservoirP.G.cells.num;
stepNPV = arroba(@NPVNet,[-1,-1,1,2],{nCells,'scale',1/10000,'sign',-1},true);

vScale = [vScale;1];
%% Instantiate the simulators for each interval, locally and for each worker.

% ss.stepClient = local (client) simulator instances 
% ss.state = scaled initial state
% ss.nv = number of algebraic variables
% ss.ci = piecewice control mapping on the client side
% ss.step =  worker simulator instances 

step = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
                                        'xScale',xScale,...
                                        'vScale',vScale,...
                                        'uScale',cellControlScales{cik},...
                                        'algFun',stepNPV,...
                                        varargin{:});
end



ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );

ss.step = step;
ss.ci = ci;



%% instantiate the objective function


%%% objective function
nCells = reservoirP.G.cells.num;
objJk = arroba(@NPVStepM,[-1,1,2,3],{nCells,'scale',1/10000,'sign',-1},true);
obj = cell(totalPredictionSteps,1);
fluid = reservoirP.fluid;
system = reservoirP.system;
for k = 1:totalPredictionSteps
    obj{k} = arroba(@mrstTimePointFuncWrapper,...
        [1,2,3],...
        {...
        objJk,...
        stepSchedules(k),...
        wellSol,...
        fluid,...
        system,...
        'xScale',...
        xScale,...
        'vScale',...
        vScale,...
        'uScale',...
        cellControlScales{callArroba(ci,{k})}...
        },true);
end



u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);


[ errorMax ] = unitTest(u{1},ss,obj,'totalSteps',3,'debug',true)


[ errorRunAdjoint ] = testRunAdjointADI( reservoirP.state, reservoirP.G, reservoirP.rock, reservoirP.fluid, reservoirP.schedule, reservoirP.system,xScale,vScale,cellControlScales)










