function ok = unitTestExample()



% Required MRST modules
mrstModule add deckformat ad-fi ad-core ad-blackoil ad-props

here = fileparts(mfilename('fullpath'));
if isempty(here)
    here = pwd();
end


% Include REMSO functionalities
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstDerivated')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'procedural')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'OOP')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'multipleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remso')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoCrossSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'singleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'utils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'testFunctions')))
addpath(genpath(fullfile(here,filesep,'reservoirData')));

%% Initialize reservoir -  the Simple reservoir
%[reservoirP] = initReservoir( 'simple10x1x10.data','Verbose',true);
[reservoirP] = initReservoir( 'reallySimpleRes.data','Verbose',true);

reservoirP.model.toleranceMB = 1e-12;
reservoirP.model.toleranceCNV = 1e-7;
reservoirP.model.toleranceWellBHP = barsa/1e-5;
reservoirP.model.toleranceWellRate = 1/day/1e-5;

model = reservoirP.model;

% do not display reservoir simulation information!
mrstVerbose off;


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
cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
    'RATE',10*meter^3/day,...
    'ORAT',10*meter^3/day,...
    'WRAT',10*meter^3/day,...
    'LRAT',10*meter^3/day,...
    'RESV',0,...
    'BHP',5*barsa));

%% instantiate the objective function as an aditional Algebraic variable


%%% The sum of the last elements in the algebraic variables is the objective
nCells = reservoirP.model.G.cells.num;
stepNPV = arroba(@NPVStepM,[1,2],{nCells,'scale',1/10000,'sign',-1},true);
%% Instantiate the simulators for each interval, locally and for each worker.

% ss.stepClient = local (client) simulator instances
% ss.state = scaled initial state
% ss.nv = number of algebraic variables
% ss.ci = piecewice control mapping on the client side
% ss.step =  worker simulator instances

step = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,stepSchedules(k),reservoirP,...
                                        'uScale',cellControlScales{cik},...
                                        'algFun',stepNPV,...
                                        varargin{:});
end
ss.state = model.toStateVector( reservoirP.state);
ss.step = step;
ss.ci = ci;


%% instantiate the objective function


%%% objective function
obj = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    obj{k} = arroba(@lastAlg,[1,2,3],{},true);
end

    u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);


[ errorMax,eCross ] = unitTest(u,ss,obj,'totalSteps',3,'debug',true);
        
ok = errorMax< 1e-3   
        
end
