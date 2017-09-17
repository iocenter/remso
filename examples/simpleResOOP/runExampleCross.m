% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global
clear classdef
clear class


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
addpath(genpath(fullfile(here,filesep,'reservoirData')));

%% Initialize reservoir -  the Simple reservoir
[reservoirP] = initReservoir( 'simple10x1x10.data','Verbose',true);
%[reservoirP] = initReservoir( 'reallySimpleRes.data','Verbose',true);

%reservoirP.model.toleranceMB = 1e-12;
%reservoirP.model.toleranceCNV = 1e-7;
%reservoirP.model.toleranceWellBHP = barsa/1e-5;
%reservoirP.model.toleranceWellRate = 1/day/1e-5;

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
targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,varargin{:});

%%  Bounds for all variables!

% Bounds for all wells!
maxProd = struct('BHP',200*barsa);
minProd = struct('BHP',(50)*barsa);
maxInj = struct('RATE',250*meter^3/day);
minInj = struct('RATE',0);


% Control input bounds for all wells!
[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProd,'minProd',minProd,...
    'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);
lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales);

% Bounds for all wells!
maxProd = struct('ORAT',200*meter^3/day,'WRAT',200*meter^3/day,'GRAT',200*meter^3/day,'BHP',200*barsa);
minProd = struct('ORAT',0,  'WRAT',0,  'GRAT',0,'BHP',(50)*barsa);
maxInj = struct('ORAT',250*meter^3/day,'WRAT',250*meter^3/day,'GRAT',250*meter^3/day,'BHP',(400)*barsa);
minInj = struct('ORAT',0,  'WRAT',0,  'GRAT',0,'BHP',(100)*barsa);

ubv = cell(totalPredictionSteps,1);
lbv = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    % during a shooting period the number of wells cannot change, sorry...
    W = stepSchedules(k).control(1).W;
    wellSol = initWellSolAD(stepSchedules(k).control(1).W,model,reservoirP.state);
    
    % wellSol bounds  (Algebraic variables bounds)
    [ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
        'maxProd',maxProd,...
        'maxInj',maxInj,...
        'minProd',minProd,...
        'minInj',minInj);
    
    ubwk = model.toWellSolVector(ubWellSol);
    lbwk = model.toWellSolVector(lbWellSol);

    
    
    
    %%%%  
    % Implement here the bound for the bounds for algFun at each shooting period
    % in this example we only have the objective (which is unbounded!)
    
    ubn = inf;
    lbn = -inf;
    
    lbv{k} = [lbwk;lbn];
    ubv{k} = [ubwk;ubn];
    
    
end

% Pay attention!, in this case the transformation to bring a mrstState to a
% stateVector is a simple scaling, that is why this works!, otherwiese it
% might be different


maxState = struct('pressure',600*barsa,'sW',1);
minState = struct('pressure',100*barsa,'sW',0.1);
ubxS = setStateValues(maxState,'nCells',nCells,'scaling',model.scaling);
lbxS = setStateValues(minState,'nCells',nCells,'scaling',model.scaling);
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Max saturation for well producer cells
for k = 1:numel(totalPredictionSteps)
    % during a shooting period the number of wells cannot change, sorry...
    W = stepSchedules(k).control(1).W;
    wellSol = initWellSolAD(W,model,reservoirP.state);
    
    prodInx  = (vertcat(wellSol.sign) < 0);
    wc    = vertcat(W(prodInx).cells);
    
    maxSat = struct('pressure',inf,'sW',1);
    ubxS = setStateValues(maxSat,'x',ubxS,'scaling',model.scaling,'cells',wc);
    ubxsatWMax = repmat({ubxS},totalPredictionSteps,1);
    ubx = cellfun(@(x1,x2)min(x1,x2),ubxsatWMax,ubx,'UniformOutput',false);
end


%% Initial Active set!
initializeActiveSet = false;
if initializeActiveSet
    vDims = cellfun(@numel,lbv);
    [ lowActive,upActive ] = activeSetFromWells(vDims,reservoirP,totalPredictionSteps);
else
    lowActive = [];
    upActive = [];
end



plotSol = @(x,u,v,xd,varargin) plotSolutionOOP( x,u,v,xd,ss,obj,model,cellControlScales,stepSchedules,controlSchedules,minState,maxState,'units','METRIC',varargin{:});

%%  Initialize from previous solution?

if exist('optimalVars.mat','file') == 2
    load('optimalVars.mat','x','u','v');
elseif exist('itVars.mat','file') == 2
    load('itVars.mat','x','u','v');
else
    x = [];
    v = [];
    u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);
    %[x] = repmat({ss.state},totalPredictionSteps,1);
end

        %% call REMSO
        
        [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
            'tol',1e-6,'lkMax',4,'debugLS',true,...
            'lowActive',lowActive,'upActive',upActive,...
            'plotFunc',plotSol,'max_iter',500,'x',x,'v',v,'debugLS',false,'saveIt',false);
        
        %% plotSolution
        plotSol(x,u,v,xd,'simFlag',true);
        