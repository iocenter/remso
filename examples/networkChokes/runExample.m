% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global

% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi ad-core ad-props

here = fileparts(mfilename('fullpath'));
if isempty(here)
here = pwd();
end

% Include REMSO functionalities
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstDerivated')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'procedural')));

addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'plottings')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'dpFunctions',filesep,'fluidProperties')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'dpFunctions',filesep,'pipeFlow')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'networkFunctions')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'auxiliaryFunctions')));

addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'multipleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remso')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoCrossSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'singleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'utils')));
addpath(genpath(fullfile(here,filesep,'reservoirData')));



%% Initialize reservoir -  the Simple reservoir
[reservoirP] = initReservoir( 'RATE10x5x10.txt','Verbose',true);

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
    W = processWells(reservoirP.G, reservoirP.rock,reservoirP.schedule.control(1),'DepthReorder', true);
end
wellSol = initWellSolLocal(W, reservoirP.state);

netSol = prodNetwork(wellSol, 'satelliteWellsNetwork', true); % instantiate the production network object

[vScale, nScale] = mrstAlg2algVar( wellSolScaling(wellSol,'bhp',5*barsa,'qWs',10*meter^3/day,'qOs',10*meter^3/day), netSolScaling(netSol, 'pressure', 5*barsa, 'flow', [], 'freq', []));

dpChokes = arroba(@chokesDp,[1,2, 3],{netSol, nScale, []}, true);

cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
    'RATE',10*meter^3/day,...
    'ORAT',10*meter^3/day,...
    'WRAT',10*meter^3/day,...
    'LRAT',10*meter^3/day,...
    'RESV',0,...
    'BHP',5*barsa));

%% instantiate the objective function as an aditional Algebraic variable


%%% The sum of the last elements in the algebraic variables is the objective
nCells = reservoirP.G.cells.num;
stepNPV = arroba(@NPVStepM,[1,2, 3],{nCells,'scale',1/10000,'sign',-1},true);
pressureScale = nScale;

vScale = [vScale; nScale;1];

[ algFun ] = concatenateMrstTargets([dpChokes, stepNPV],false, [1; 1]);

step = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
                                        'xScale',xScale,...
                                        'vScale',vScale,...
                                        'uScale',cellControlScales{cik},...
                                        'algFun',algFun,...
                                        varargin{:});
end

ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
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
minProd = struct('BHP',100*barsa);
maxProd = struct('BHP',500*barsa);
minInj = struct('RATE',5*meter^3/day);
maxInj = struct('RATE',300*meter^3/day);


% Control input bounds for all wells!
[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProd,'minProd',minProd,...
    'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);

lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales);

% Bounds for all wells!
minProd = struct('ORAT', 5*meter^3/day,  'WRAT', -inf,  'GRAT', -inf,'BHP',-inf);
maxProd = struct('ORAT', inf,'WRAT', inf,'GRAT', inf,'BHP',inf);
minInj = struct('ORAT',-inf,  'WRAT', -inf,  'GRAT', -inf,'BHP', -inf);
maxInj = struct('ORAT',inf,'WRAT', inf ,'GRAT', inf,'BHP',800*barsa);

% wellSol bounds  (Algebraic variables bounds)
[ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
    'maxProd',maxProd,...
    'maxInj',maxInj,...
    'minProd',minProd,...
    'minInj',minInj);
ubvS = wellSol2algVar(ubWellSol,'vScale',vScale);
lbvS = wellSol2algVar(lbWellSol,'vScale',vScale);
lbv = repmat({[lbvS;  0.*nScale; -inf]},totalPredictionSteps,1);
ubv = repmat({[ubvS;   inf.*nScale; inf]},totalPredictionSteps,1);

% State lower and upper - bounds
maxState = struct('pressure',1000*barsa,'sW',1);
minState = struct('pressure',50*barsa,'sW',0.1);

ubxS = setStateValues(maxState,'nCells',nCells,'xScale',xScale);
lbxS = setStateValues(minState,'nCells',nCells,'xScale',xScale);
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Max saturation for well producer cells

% get well perforation grid block set
prodInx  = (vertcat(wellSol.sign) < 0);
wc    = vertcat(W(prodInx).cells);

maxSat = struct('pressure',inf,'sW',1);
%ubxS = setStateValues(maxSat,'x',ubxS,'xScale',xScale,'cells',wc,'nCells',nCells);
ubxS = setStateValues(maxSat,'x',ubxS,'xScale',xScale,'cells',wc);
ubxsatWMax = repmat({ubxS},totalPredictionSteps,1);
ubx = cellfun(@(x1,x2)min(x1,x2),ubxsatWMax,ubx,'UniformOutput',false);

%% A plot function to display information at each iteration

times.steps = [stepSchedules(1).time;arrayfun(@(x)(x.time+sum(x.step.val))/day,stepSchedules)];
times.tPieceSteps = cell2mat(arrayfun(@(x)[x;x],times.steps,'UniformOutput',false));
times.tPieceSteps = times.tPieceSteps(2:end-1);

times.controls = [controlSchedules(1).time;arrayfun(@(x)(x.time+sum(x.step.val))/day,controlSchedules)];
times.tPieceControls = cell2mat(arrayfun(@(x)[x;x],times.controls,'UniformOutput',false));
times.tPieceControls = times.tPieceControls(2:end-1);

cellControlScalesPlot = schedules2CellControls(schedulesScaling( controlSchedules,'RATE',1/(meter^3/day),...
    'ORAT',1/(meter^3/day),...
    'WRAT',1/(meter^3/day),...
    'LRAT',1/(meter^3/day),...
    'RESV',0,...
    'BHP',1/barsa));

[uMlb] = scaleSchedulePlot(lbu,controlSchedules,cellControlScales,cellControlScalesPlot);
[uLimLb] = min(uMlb,[],2);
ulbPlob = cell2mat(arrayfun(@(x)[x,x],uMlb,'UniformOutput',false));


[uMub] = scaleSchedulePlot(ubu,controlSchedules,cellControlScales,cellControlScalesPlot);
[uLimUb] = max(uMub,[],2);
uubPlot = cell2mat(arrayfun(@(x)[x,x],uMub,'UniformOutput',false));


% be carefull, plotting the result of a forward simulation at each
% iteration may be very expensive!
% use simFlag to do it when you need it!
simFunc =@(sch,varargin) runScheduleADI(reservoirP.state, reservoirP.G, reservoirP.rock, reservoirP.system, sch,'force_step',false,varargin{:});


wc    = vertcat(W.cells);
fPlot = @(x)[max(x);min(x);x(wc)];

plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d, lbv, ubv, lbu, ubu, ss,obj,times,xScale,cellControlScales,vScale, nScale, ...
    cellControlScalesPlot,controlSchedules,wellSol, netSol, ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true, 'plotNetsol', true, ...
    'numNetConstraints', numel(nScale), 'plotNetControls', false, 'pressureCst', numel(pressureScale), 'plotChokeConstraints', true, ...
    'plotSchedules',false,'pF',fPlot,'sF',fPlot, 'plotCumulativeObjective', true, 'reservoirP', reservoirP, 'plotNetwork', false, 'dpFunction', @dpBeggsBrillJDJ,  varargin{:});

%%  Initialize from previous solution?
loadPrevSolution = true;
optimize = true;
plotSolution = false;

if loadPrevSolution
    load('itVars.mat','x','u','v');
else
    x = [];
    v = [];
    u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);
end
   
%% call REMSO
if optimize
    [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
                    'tol',1e-6,'lkMax',4,'debugLS',true,...
                    'lowActive',[],'upActive',[],...
                    'plotFunc',plotSol,'max_iter', 500,'x',x,'v',v,'debugLS',false,'saveIt',true, 'computeCrossTerm', false, 'condense', true);
end

%% plotSolution
if plotSolution
    if ~optimize && ~loadPrevSolution
        [~, ~, ~, simVars, x, v] = simulateSystemSS(u, ss, []);
    end
    xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
    plotSol(x,u,v,xd, 'simFlag', false)     
end
      