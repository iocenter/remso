% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global
clear classdef
clear class


runInParallel = false;

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
% addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'OOP')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'utils')));

addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'multipleShooting')));
if runInParallel
    addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'parallel')));    
end
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remso')));

if ~runInParallel
    addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoSequential')));        
end

addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoCrossSequential')));

if runInParallel
    addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoCross')));
end



addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'singleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'utils')));
addpath(genpath(fullfile(here,filesep,'reservoirData')));


addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'plottings')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'networkFunctions')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'auxiliaryFunctions')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'dpFunctions',filesep,'fluidProperties')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'dpFunctions',filesep,'pipeFlow')));

% Open a matlab pool depending on the machine availability
if runInParallel
    initPool('restart',true);
end


%% Initialize reservoir
% [reservoirP] = initReservoir('RATE10x5x10.txt', 'Verbose',true);
[reservoirP] = loadEgg('./reservoirData/Egg_Model_Data_Files_v2/MRST');

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
[vScale] = mrstAlg2algVar( wellSolScaling(wellSol,'bhp',5*barsa,'qWs',10*meter^3/day,'qOs',10*meter^3/day, 'freq', 15), netSolScaling(netSol));

freqScale = [15;15;15;15]; % in Hz
flowScale = [5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day); ...
    5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day)];

pressureScale = [5*barsa;5*barsa;5*barsa;5*barsa;];

% number of pump stages
numStages =  [90; 90; 90; 90];
% bounds for flowing rates through the pump at 60 Hz

%% Infeasible initial guess
% qlMin = [105*(meter^3/day); ... % PROD1
%          65*(meter^3/day); ...% PROD3
%          225*(meter^3/day); ... % PROD2
%          95*(meter^3/day)];    % PROD4
%      
% qlMax = [215*(meter^3/day); ...  % PROD1
%          195*(meter^3/day); ...  % PROD3
%          405*(meter^3/day); ...  % PROD2
%          205*(meter^3/day);];    % PROD4

%% Feasible initial guess
qlMin = [85*(meter^3/day); ... % PROD1
         65*(meter^3/day); ...% PROD3
         225*(meter^3/day); ... % PROD2
         75*(meter^3/day)];    % PROD4
     
qlMax = [195*(meter^3/day); ...  % PROD1
         195*(meter^3/day); ...  % PROD3
         405*(meter^3/day); ...  % PROD2
         185*(meter^3/day);];    % PROD4\

% bounds for pump frequencies in Hz
freqMin = [40; ... % PROD1
           40; ... % PROD3
           40; ... % PROD2
           40;];   % PROD4
       
freqMax = [80; ... % PROD1
           80; ... % PROD3
           80; ... % PROD2
           80;];   % PROD4
       
baseFreq = [60; 60; 60; 60;]; % in Hz


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

networkJointObj = arroba(@networkJointNPVConstraints,[1,2],{nCells, netSol, freqScale, pressureScale, flowScale, numStages, baseFreq, qlMin, qlMax, 'scale',1/100000,'sign',-1, 'dpFunction', @dpBeggsBrillJDJ, 'finiteDiff', true, 'forwardGradient', true, 'extremePoints', extremePoints},true);

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
minProd = struct('BHP',100*barsa, 'ORAT',  1*meter^3/day);

% maxProd = struct('BHP',200*barsa, 'ORAT', 220*meter^3/day); original val
maxProd = struct('BHP',450*barsa, 'ORAT', inf*meter^3/day);

% minInj = struct('RATE',100*meter^3/day); % original val
minInj = struct('RATE',1*meter^3/day);
% maxInj = struct('RATE',300*meter^3/day); original val

% maxInj = struct('RATE',300*meter^3/day);
maxInj = struct('RATE',500*meter^3/day);

% Control input bounds for all wells!

[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProd,'minProd',minProd,...
    'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);
lbw = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);
ubw = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);


cellControlScale = cellfun(@(wi) wi,cellControlScales,'uniformOutput', false);

lbu = cellfun(@(wi) wi,lbw, 'UniformOutput',false);
ubu = cellfun(@(wi) wi,ubw, 'UniformOutput',false);


% Bounds for all wells!
% minProd = struct('ORAT',1*meter^3/day,  'WRAT',1*meter^3/day,  'GRAT',
% -inf,'BHP',130*barsa); original val
minProd = struct('ORAT', -inf,  'WRAT', -inf,  'GRAT', -inf,'BHP',-inf);
% maxProd = struct('ORAT',220*meter^3/day,'WRAT',150*meter^3/day,'GRAT',
% inf,'BHP',350*barsa); original val
maxProd = struct('ORAT', inf,'WRAT', inf,'GRAT', inf,'BHP',inf);

% minInj = struct('ORAT',-inf,  'WRAT',100*meter^3/day,  'GRAT',
% -inf,'BHP', 5*barsa); original val
minInj = struct('ORAT',-inf,  'WRAT', -inf,  'GRAT', -inf,'BHP', -inf);
% maxInj = struct('ORAT',inf,'WRAT',300*meter^3/day,'GRAT',
% inf,'BHP',500*barsa); original val
maxInj = struct('ORAT',inf,'WRAT', inf ,'GRAT', inf,'BHP', inf);

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

% Non-Linear Pump Map Constraints
lbv = repmat({[lbvS; 0./flowScale;   freqMin./freqScale;   -inf]},totalPredictionSteps,1);
ubv = repmat({[ubvS; inf./flowScale; freqMax./freqScale;    inf]},totalPredictionSteps,1);

% lbv = repmat({[lbvS; -inf./flowScale;   -inf./freqScale;   -inf]},totalPredictionSteps,1);
% ubv = repmat({[ubvS; inf./flowScale;    inf./freqScale;    inf]},totalPredictionSteps,1);

% State lower and upper - bounds
maxState = struct('pressure',1000*barsa,'sW',1);
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

cellControlScalesPlot = cellfun(@(w) w, cellControlScalesPlot, 'UniformOutput',false);

cellControlScales  = cellfun(@(w) w , cellControlScales ,'uniformOutput', false);

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
    cellControlScalesPlot,controlSchedules,wellSol, netSol, ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true, 'plotNetsol', false, ...
    'numNetConstraints', numel(nScale), 'plotNetControls', false, 'freqCst', numel(freqScale), 'pressureCst',numel(pressureScale),  'flowCst',numel(flowScale), ...
    'plotSchedules',false,'pF',fPlot,'sF',fPlot, 'fixedWells', fixedWells, 'extremePoints', extremePoints, 'plotCumulativeObjective', true, 'qlMin', qlMin,  'qlMax', qlMax, 'nStages', numStages, ...
    'freqMin', freqMin, 'freqMax', freqMax, 'baseFreq', baseFreq, 'reservoirP', reservoirP, 'plotNetwork', true, 'wc', true, 'dpFunction', @dpBeggsBrillJDJ, varargin{:});


% remove network control to initialize well controls vector (w)
cellControlScales = cellfun(@(w) w(1:end) ,cellControlScales, 'UniformOutput', false);

%%  Initialize from previous solution?

x = [];
v = [];
w  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);

cellControlScales = cellfun(@(w) w , cellControlScales ,'uniformOutput', false);
u = cellfun(@(wi)wi, w,'UniformOutput',false);
cellControlScalesPlot = cellfun(@(w) w, cellControlScalesPlot,'uniformOutput', false);

controlWriter = @(u,i) controlWriterMRST(u,i,controlSchedules,cellControlScales,'filename',['./controls/schedule' num2str(i) '.inc'], 'fixedWells', fixedWells);

loadPrevSolution = false;
optimize = false;
loadSingleControl = false;
load10Controls = true;
plotSolution = true;
saveFigures = true;

if loadPrevSolution
   load itVars;
end

if loadSingleControl
   load singleControlOpt;    
   u = repmat(u, numel(ubu),1); %% refine u    
end

if load10Controls
   load 10controlsConstrained; 
end

if optimize
    warning('OFF','ALL');    
    tic
    [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
        'skipRelaxRatio',inf,'tol',1e-4,'lkMax',4, ...
        'lowActive',lowActive,'upActive',upActive,...
        'plotFunc',plotSol,'max_iter', 500,'x',x,'v',v,'debugLS',false,'saveIt',true, 'computeCrossTerm', false, 'condense', true,'controlWriter',controlWriter, 'qpAlgorithm', 2);
    compTime = toc;
    save('time.mat', 'compTime');
end

if  plotSolution
    if optimize
        plotSol(x,u,v,xd, 'simFlag', false);
    elseif loadPrevSolution || loadSingleControl || load10Controls
        xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
        plotSol(x,u,v,xd, 'simFlag', false)
    else
        [~, ~, ~, simVars, x, v] = simulateSystemSS(u, ss, [])
        xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
        plotSol(x,u,v,xd, 'simFlag', false);
    end
    
    if saveFigures
        figlist=findobj('type','figure');
        dirname = 'figs/';
        if ~(exist(dirname,'dir')==7)
            mkdir(dirname);
        end
        
        for i=1:numel(figlist)
            saveas(figlist(i),fullfile(dirname, ['figure' num2str(figlist(i).Number) '.eps']), 'epsc');
        end
        close all;
    end
end
