% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global

runInParallel = false;
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

% Open a matlab pool depending on the machine availability
if runInParallel
initPool('restart',true);
end

%% Initialize reservoir -  the Simple reservoir
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

freqScale = []; % in Hz
flowScale = [];
pressureScale = [];

%% network controls
pScale = [];
p  = [];

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

stepNPV = arroba(@NPVStepM,[1,2],{nCells, ...
                  'scale',1e-07, 'OilPrice', 100, 'WaterProductionCost', 5, ...
                   'WaterInjectionCost', 0.1, 'DiscountFactor', 0, 'sign',-1, ...
                  },true);

vScale = [vScale; 1];

stepClient = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    stepClient{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
        'xScale',xScale,...
        'vScale',vScale,...
        'uScale',cellControlScales{cik},...
        'algFun',stepNPV,...
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
            'algFun',stepNPV...
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

maxInj = struct('RATE',500*meter^3/day);

% Control input bounds for all wells!

[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProd,'minProd',minProd,...
    'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);
lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);


% Bounds for all wells!
minProd = struct('ORAT', -inf*meter^3/day,  'WRAT', -inf,  'GRAT', -inf,'BHP',-inf);
maxProd = struct('ORAT', inf,'WRAT', inf,'GRAT', inf,'BHP',inf);

minInj = struct('ORAT',-inf,  'WRAT', -inf,  'GRAT', -inf,'BHP', -inf);
maxInj = struct('ORAT',inf,'WRAT', inf ,'GRAT', inf,'BHP',inf);

% wellSol bounds  (Algebraic variables bounds)
[ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
    'maxProd',maxProd,...
    'maxInj',maxInj,...
    'minProd',minProd,...
    'minInj',minInj);

ubvS = wellSol2algVar(ubWellSol,'vScale',vScale);
lbvS = wellSol2algVar(lbWellSol,'vScale',vScale);

%% without network constraints
lbv = repmat({[lbvS;  -inf]},totalPredictionSteps,1);
ubv = repmat({[ubvS;  inf]},totalPredictionSteps,1);

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

plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d, lbv, ubv, lbu, ubu, ss,objClient,times,xScale,cellControlScales,vScale, [], ...
    cellControlScalesPlot,controlSchedules,wellSol, netSol, ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true, 'plotNetsol', false, ...
    'numNetConstraints', [], 'plotNetControls', false, 'numNetControls', numel(pScale), 'freqCst', numel(freqScale), 'pressureCst',numel(pressureScale),  'flowCst',numel(flowScale), ...
    'plotSchedules',false,'pF',fPlot,'sF',fPlot, 'fixedWells', fixedWells, 'extremePoints', [], 'plotCumulativeObjective', true, 'qlMin', [],  'qlMax', [], 'nStages', [], ...
    'freqMin', [], 'freqMax', [], 'baseFreq', [], 'reservoirP', reservoirP, 'plotNetwork', false, 'wc', true, 'dpFunction', @dpBeggsBrillJDJ,  varargin{:});

% remove network control to initialize well controls vector (w)
cellControlScales = cellfun(@(w) w(1:end-numel(p)) ,cellControlScales, 'UniformOutput', false);

%%  Initialize from previous solution?
x = [];
v = [];
u = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);

controlWriter = @(u,i) controlWriterMRST(u,i,controlSchedules,cellControlScales,'filename',['./controls/schedule' num2str(i) '.inc'], 'fixedWells', fixedWells);

loadPrevSolution = false;
optimize = true;
plotSolution = false;
        if loadPrevSolution
            load itVars;
        end
        
        if optimize
            tic
            [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
                'tol',1e-4,'lkMax',4,'debugLS',true,...
				'skipRelaxRatio',inf,...
                'lowActive',lowActive,'upActive',upActive,...
                'plotFunc',plotSol,'max_iter', 500,'x',x,'v',v,'debugLS',false,'saveIt',true, 'computeCrossTerm', false, 'condense', true,'controlWriter',controlWriter, 'qpAlgorithm', 2);
            compTime = toc;
            save('time.mat', 'compTime');
        end
        
        if  plotSolution
            if optimize
                plotSol(x,u,v,xd, 'simFlag', false);
            elseif loadPrevSolution
                xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
                plotSol(x,u,v,xd, 'simFlag', false)
            else
                [~, ~, ~, simVars, x, v] = simulateSystemSS(u, ss, [])
                xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
                plotSol(x,u,v,xd, 'simFlag', false)
            end
        end        