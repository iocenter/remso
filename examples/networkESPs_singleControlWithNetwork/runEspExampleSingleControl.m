% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global


runInParallel = true;

% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi ad-core ad-props

% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/procedural'));

addpath(genpath('../../netLink'));
addpath(genpath('../../netLink/plottings'));
addpath(genpath('../../netLink/dpFunctions/fluidProperties'));
addpath(genpath('../../netLink/dpFunctions/pipeFlow'));
addpath(genpath('../../netLink/networkFunctions'));
addpath(genpath('../../netLink/auxiliaryFunctions'));

addpath(genpath('../../optimization/multipleShooting'));
if runInParallel
    addpath(genpath('../../optimization/parallel'));
end
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
if ~runInParallel
addpath(genpath('../../optimization/remsoSequential'));
end
addpath(genpath('../../optimization/remsoCrossSequential'));
if runInParallel
    addpath(genpath('../../optimization/remsoCross'));
end
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));

% Open a matlab pool depending on the machine availability
if runInParallel
initPool('restart',true);
end


%% Initialize reservoir -  the Simple reservoir
[reservoirP] = initReservoir('RATE10x5x10.txt', 'Verbose',true);


reservoirP.schedule.control = reservoirP.schedule.control(1);
reservoirP.schedule.step.control = ones(size(reservoirP.schedule.step.control));


% do not display reservoir simulation information!
mrstVerbose off;

% Number of reservoir grid-blocks
nCells = reservoirP.G.cells.num;

%% Multiple shooting problem set up
totalPredictionSteps = numel(reservoirP.schedule.step.val);  % MS intervals

% Schedule partition for each control period and for each simulated step
lastControlSteps = findControlFinalSteps( reservoirP.schedule.step.control );
controlSchedules = multipleSchedules(reservoirP.schedule,lastControlSteps);

uUnscaled  = schedules2CellControls( controlSchedules);
uDims = cellfun(@(uu)size(uu,1),uUnscaled);
totalControlSteps = length(uUnscaled);

stepSchedules = multipleSchedules(reservoirP.schedule,1:totalPredictionSteps);

% Piecewise linear control -- mapping the step index to the corresponding
% control
ci  = arroba(@controlIncidence, 2 ,{reservoirP.schedule.step.control});


if runInParallel
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


[workerCondensingSchedule,clientCondensingSchedule,uStart,workerLoad,avgW] = divideCondensingLoad(totalPredictionSteps,ci,uDims,nWorkers);


jobSchedule.clientCondensingSchedule = clientCondensingSchedule;
jobSchedule.workerCondensingSchedule = workerCondensingSchedule;
jobSchedule.uStart = uStart;
end


%% Variables Scaling
xScale = setStateValues(struct('pressure',5*barsa,'sW',0.01),'nCells',nCells);

W =  reservoirP.schedule.control.W;
wellSol = initWellSolLocal(W, reservoirP.state);
for k = 1:numel(wellSol)
    wellSol(k).qGs = 0;
end
nW = numel(W);

%% Fixing Injectors
%     fixedWells = find(vertcat(W.sign) == 1); % fixing injectors
fixedWells = [];
controlWells = setdiff(1:nW, fixedWells);

% Instantiate the production network object
netSol = prodNetwork(wellSol, 'espNetwork', true, 'withPumps', true);

%%TODO: separate scalling of vk and nk.
%% Scallings
[vScale, freqScale] = mrstAlg2algVar( wellSolScaling(wellSol,'bhp',5*barsa,'qWs',10*meter^3/day,'qOs',10*meter^3/day, 'freq', 15), netSolScaling(netSol));

%     freqScale = [];
%     flowScale = [];
freqScale = [15;15;15;15;15]; % in Hz
flowScale = [5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day);5*(meter^3/day); ...
    5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day);5*(meter^3/day)];

pressureScale = [5*barsa;5*barsa;5*barsa;5*barsa;5*barsa];

%% network controls
pScale = [];
p  = [];

% number of pump stages
numStages =  [70; 70; 70; 70; 70];
% bounds for flowing rates through the pump at 60 Hz
qlMin = [45*(meter^3/day); 45*(meter^3/day); 45*(meter^3/day); 45*(meter^3/day); 45*(meter^3/day)];
qlMax = [100*(meter^3/day); 100*(meter^3/day); 100*(meter^3/day); 100*(meter^3/day); 100*(meter^3/day)];

% bounds for pump frequencies
freqMin = [30; 30; 30; 30; 30]; % in Hz
freqMax = [90; 90; 90; 90; 90]; % in Hz
baseFreq = [60; 60; 60; 60; 60]; % in Hz


%% run network four times to obtain the extreme points of the pumps maps
qminFmin = pump_rate(freqMin, qlMin, baseFreq);
qminFmax = pump_rate(freqMax, qlMin, baseFreq);
qmaxFmin = pump_rate(freqMin, qlMax, baseFreq);
qmaxFmax = pump_rate(freqMax, qlMax, baseFreq);

qf = cell(4,1);
dp = cell(4,1);

qf{1}  = qminFmin./(meter^3/day); qf{2} = qminFmax./(meter^3/day);
qf{3} = qmaxFmin./(meter^3/day); qf{4} = qmaxFmax./(meter^3/day);

str = netSol.E(1).stream; % network has default stream for subsea pipeline

oilDens = str.oil_dens;
highestDens = 0.6*str.water_dens + 0.4*str.oil_dens; % maximum of 0.6 water cut (based on experiments, not limited in practice)

dp{1} = calcDp(qminFmin, freqMin, baseFreq, numStages, 'mixDensity',  oilDens);
dp{2} = calcDp(qminFmax, freqMax, baseFreq, numStages, 'mixDensity', highestDens);
dp{3} = calcDp(qmaxFmin, freqMin, baseFreq, numStages, 'mixDensity', oilDens);
dp{4} = calcDp(qmaxFmax, freqMax, baseFreq, numStages, 'mixDensity', highestDens);

extremePoints = cell(4,1);
for i=1:4
    extremePoints{i} = [qf{i}, dp{i}];
end
extremePoints = [];

cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
    'RATE',10*meter^3/day,...
    'ORAT',10*meter^3/day,...
    'WRAT',10*meter^3/day,...
    'LRAT',10*meter^3/day,...
    'RESV',0,...
    'BHP',5*barsa),'fixedWells', fixedWells);

%% instantiate the objective function as an aditional Algebraic variable
%%% The sum of the last elements in the algebraic variables is the objective
nCells = reservoirP.G.cells.num;
nScale  = [flowScale; freqScale; pressureScale];
vScale = [vScale; nScale; 1];

networkJointObj = arroba(@networkJointNPVConstraints,[1,2, 3],{nCells, netSol, freqScale, pressureScale, flowScale, numStages, baseFreq, qlMin, qlMax, pScale,   'scale',1/100000,'sign',-1, 'dpFunction', @dpBeggsBrillJDJ, 'finiteDiff', true, 'forwardGradient', true, 'extremePoints', extremePoints},true);

%% Instantiate the simulators for each interval, locally and for each worker.
stepClient = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    stepClient{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
        'xScale',xScale,...
        'vScale',vScale,...
        'uScale',cellControlScales{cik},...
        'algFun',networkJointObj,...
        'fixedWells', fixedWells, ...
        'saveTargetJac', true,...
        varargin{:});
    
end

if runInParallel
spmd
    
    stepW = cell(totalPredictionSteps,1);
    for k=1:totalPredictionSteps
        cik = callArroba(ci,{k});
        stepW{k} = arroba(@mrstStep,...
            [1,2],...
            {...
            @mrstSimulationStep,...
            wellSol,...
            stepSchedules(k),...
            reservoirP,...
            'xScale',...
            xScale,...
            'vScale',...
            vScale,...
            'uScale',...
            cellControlScales{cik},...
            'algFun',networkJointObj...
            'fixedWells', fixedWells, ...
            'saveTargetJac', true,...
            },...
            true);
    end
    stepC = stepW;
end
end
ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
if runInParallel
ss.jobSchedule = jobSchedule;
ss.work2Job = work2Job;
ss.step = stepC;
ss.stepClient = stepClient;
else
ss.step = stepClient;
end
ss.ci = ci;

%% instantiate the objective function
%%% objective function on the client side (for plotting!)
objClient = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    objClient{k} = arroba(@lastAlg,[1,2,3],{},true);
end
%%% Investigate if it is efficient to evaluate the objective in the workers
if runInParallel
spmd
    nJobsW = numel(work2Job);
    objW = cell(nJobsW,1);
    for i = 1:nJobsW
        objW{i} = arroba(@lastAlg,[1,2,3],{},true);
    end
    obj = objW;
end
targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,jobSchedule,work2Job,varargin{:});
else
obj = objClient;
targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,varargin{:});
end

%%  Bounds for all variables!

% Bounds for all wells!
% minProd = struct('BHP',130*barsa, 'ORAT', 1*meter^3/day); original val
minProd = struct('BHP',100*barsa, 'ORAT',  5*meter^3/day, 'FREQ', 15);

% maxProd = struct('BHP',200*barsa, 'ORAT', 220*meter^3/day); original val
maxProd = struct('BHP',500*barsa, 'ORAT', 200*meter^3/day, 'FREQ', 100);

% minInj = struct('RATE',100*meter^3/day); % original val
minInj = struct('RATE',5*meter^3/day);
% maxInj = struct('RATE',300*meter^3/day); original val

maxInj = struct('RATE',300*meter^3/day);

% Control input bounds for all wells!

[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProd,'minProd',minProd,...
    'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);
lbw = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);
ubw = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);


cellControlScale = cellfun(@(wi) [wi; pScale],cellControlScales,'uniformOutput', false);

lbu = cellfun(@(wi)[wi; 5*barsa./pScale],lbw, 'UniformOutput',false);
ubu = cellfun(@(wi)[wi; 30*barsa./pScale],ubw, 'UniformOutput',false);


% Bounds for all wells!
% minProd = struct('ORAT',1*meter^3/day,  'WRAT',1*meter^3/day,  'GRAT',
% -inf,'BHP',130*barsa); original val
minProd = struct('ORAT', 5*meter^3/day,  'WRAT', 0*meter^3/day,  'GRAT', -inf,'BHP',100*barsa);
% maxProd = struct('ORAT',220*meter^3/day,'WRAT',150*meter^3/day,'GRAT',
% inf,'BHP',350*barsa); original val
maxProd = struct('ORAT',200*meter^3/day,'WRAT',200*meter^3/day,'GRAT', inf,'BHP',500*barsa);

% minInj = struct('ORAT',-inf,  'WRAT',100*meter^3/day,  'GRAT',
% -inf,'BHP', 5*barsa); original val
minInj = struct('ORAT',-inf,  'WRAT', 5*meter^3/day,  'GRAT', -inf,'BHP', 5*barsa);
% maxInj = struct('ORAT',inf,'WRAT',300*meter^3/day,'GRAT',
% inf,'BHP',500*barsa); original val
maxInj = struct('ORAT',inf,'WRAT', 300*meter^3/day,'GRAT', inf,'BHP',800*barsa);

% wellSol bounds  (Algebraic variables bounds)
[ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
    'maxProd',maxProd,...
    'maxInj',maxInj,...
    'minProd',minProd,...
    'minInj',minInj);

ubvS = wellSol2algVar(ubWellSol,'vScale',vScale);
lbvS = wellSol2algVar(lbWellSol,'vScale',vScale);

%% Linear Approx. of Pump Map
% lbv = repmat({[lbvS; 0./flowScale;   0*barsa./pressureScale; 0*barsa./pressureScale;  -inf*barsa./pressureScale; -inf]},totalPredictionSteps,1);
% ubv = repmat({[ubvS; inf./flowScale; inf*barsa./pressureScale; inf*barsa./pressureScale;  0*barsa./pressureScale; inf]},totalPredictionSteps,1);

%% Non-Linear Pump Map Constraints
lbv = repmat({[lbvS; 0./flowScale;   freqMin./freqScale; -inf]},totalPredictionSteps,1);
ubv = repmat({[ubvS; inf./flowScale; freqMax./freqScale;  inf]},totalPredictionSteps,1);

% State lower and upper - bounds
maxState = struct('pressure',500*barsa,'sW',1);
minState = struct('pressure',50*barsa,'sW',0.1);
ubxS = setStateValues(maxState,'nCells',nCells,'xScale',xScale);
lbxS = setStateValues(minState,'nCells',nCells,'xScale',xScale);
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lowActive = [];
upActive = [];

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

cellControlScalesPlot = cellfun(@(w) [w;pScale], cellControlScalesPlot, 'UniformOutput',false);

cellControlScales  = cellfun(@(w) [w; pScale] , cellControlScales ,'uniformOutput', false);

[uMlb] = scaleSchedulePlot(lbu,controlSchedules,cellControlScales,cellControlScalesPlot, 'fixedWells', fixedWells);
[uLimLb] = min(uMlb,[],2);
ulbPlob = cell2mat(arrayfun(@(x)[x,x],uMlb,'UniformOutput',false));


[uMub] = scaleSchedulePlot(ubu,controlSchedules,cellControlScales,cellControlScalesPlot, 'fixedWells', fixedWells);
[uLimUb] = max(uMub,[],2);
uubPlot = cell2mat(arrayfun(@(x)[x,x],uMub,'UniformOutput',false));


% be carefull, plotting the result of a forward simulation at each
% iteration may be very expensive!
% use simFlag to do it when you need it!
simFunc =@(sch,varargin) runScheduleADI(reservoirP.state, reservoirP.G, reservoirP.rock, reservoirP.system, sch,'force_step',false,varargin{:});

wc    = vertcat(W.cells);
fPlot = @(x)[max(x);min(x);x(wc)];

plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d, lbv, ubv, lbu, ubu, ss,objClient,times,xScale,cellControlScales,vScale, nScale, ...
    cellControlScalesPlot,controlSchedules,wellSol, netSol, ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true, 'plotNetsol', true, ...
    'numNetConstraints', numel(nScale), 'plotNetControls', false, 'numNetControls', numel(pScale), 'freqCst', numel(freqScale), 'pressureCst',numel(pressureScale),  'flowCst',numel(flowScale), ...
    'plotSchedules',false,'pF',fPlot,'sF',fPlot, 'fixedWells', fixedWells, 'extremePoints', extremePoints, 'plotCumulativeObjective', true, 'qlMin', qlMin,  'qlMax', qlMax, 'nStages', numStages, ...
    'freqMin', freqMin, 'freqMax', freqMax, 'baseFreq', baseFreq, 'reservoirP', reservoirP, 'plotNetwork', true, 'wc', true, 'dpFunction', @dpBeggsBrillJDJ, varargin{:});


% remove network control to initialize well controls vector (w)
cellControlScales = cellfun(@(w) w(1:end-numel(p)) ,cellControlScales, 'UniformOutput', false);

%%  Initialize from previous solution?

x = [];
v = [];
w  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);

cellControlScales = cellfun(@(w) [w; pScale] , cellControlScales ,'uniformOutput', false);
u = cellfun(@(wi)[wi;p],w,'UniformOutput',false);

cellControlScalesPlot = cellfun(@(w) [w; pScale], cellControlScalesPlot,'uniformOutput', false);

controlWriter = @(u,i) controlWriterMRST(u,i,controlSchedules,cellControlScales,'filename',['./controls/schedule' num2str(i) '.inc'], 'fixedWells', fixedWells);

loadPrevSolution = false;
optimize = true;
plotSolution = false;

if loadPrevSolution
    load itVars;
end

if optimize
    [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
        'tol',1e-4,'lkMax',4,'debugLS',true,...
        'skipRelaxRatio',inf,...
        'lowActive',lowActive,'upActive',upActive,...
        'plotFunc',plotSol,'max_iter', 500,'x',x,'v',v,'debugLS',false,'saveIt',true, 'computeCrossTerm', false, 'condense', true,'controlWriter',controlWriter);
end

if  plotSolution
    if optimize
        plotSol(x,u,v,xd, 'simFlag', false);
    elseif ~loadPrevSolution
        xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
        plotSol(x,u,v,xd, 'simFlag', false)
    else
        [~, ~, ~, simVars, x, v] = simulateSystemSS(u, ss, [])
        xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
        plotSol(x,u,v,xd, 'simFlag', false)
    end
end
