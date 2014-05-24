% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

%{

This script instantiate an optimal control problem and solve it with REMSO.

The problem is based on the Egg model. Please donwnload the Egg Model
instance from:

http://dx.doi.org/10.4121/uuid:916c86cd-3558-4672-829a-105c62985ab2

and place the files related to MRST in:

./reservoirData/Egg_Model_Data_Files_v2/MRST


%}



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
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));

% Open a matlab pool depending on the machine availability
initPool();


%% Initialize reservoir the Egg model
[reservoirP] = loadEgg('./reservoirData/Egg_Model_Data_Files_v2/MRST');

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
ci = @(k) controlIncidence(reservoirP.schedule.step.control,k);
ciW  = arroba(@controlIncidence,2,{reservoirP.schedule.step.control});


%%  Who will do what - Distribute the computational effort!
nWorkers = getNumWorkers() ;
if nWorkers == 0
    nWorkers = 1;
end
[ jobSchedule ] = divideJobsSequentially(totalPredictionSteps ,nWorkers);
jobSchedule.nW = nWorkers;

work2Job = Composite();
for w = 1:nWorkers
    work2Job{w} = jobSchedule.work2Job{w};
end


%% Variables Scaling
xScale = stateMrst2stateVector( stateScaling(reservoirP.state,...
    'pressure',5*barsa,...
    's',0.01) );


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


%% Instantiate the simulators for each interval, locally and for each worker.

% ss.stepClient = local (client) simulator instances 
% ss.state = scaled initial state
% ss.nv = number of algebraic variables
% ss.ci = piecewice control mapping on the client side
% ss.ciW = piecewice control mapping on the workers side
% ss.jobSchedule = Step distribution among workers client side
% ss.work2Job = Step distribution among workers worker side
% ss.step =  worker simulator instances 

stepClient = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    ss.stepClient{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,'xScale',xScale,'vScale',vScale,'uScale',cellControlScales{cik},varargin{:});
end


stepSchedulesW = distributeVariables(stepSchedules,jobSchedule);
spmd
    
    nJobsW = numel(work2Job);
    stepW = cell(nJobsW,1);
    for i = 1:nJobsW
        k = work2Job(i);
        cik = callArroba(ciW,{k});
        stepW{i} = arroba(@mrstStep,...
            [1,2],...
            {...
            @mrstSimulationStep,...
            wellSol,...
            stepSchedulesW(i),...
            reservoirP,...
            'xScale',...
            xScale,...
            'vScale',...
            vScale,...
            'uScale',...
            cellControlScales{cik}...
            },...
            true);
    end
    step = stepW;
end

ss.state = stateMrst2stateVector( reservoirP.state,'doScale',true,'xScale',xScale );
ss.nv = numel(vScale);
ss.jobSchedule = jobSchedule;
ss.work2Job = work2Job;
ss.step = step;
ss.ci = ci;
ss.ciW = ciW;



%% instantiate the objective function


% the objective function is a separable, exploit this!
nCells = reservoirP.G.cells.num;
objJk = arroba(@NPVStepM,[-1,1,2],{nCells,'scale',1/100000,'sign',-1,'WaterProductionCost',0.01},true);
fW = @mrstTimePointFuncWrapper;
spmd
    objW = cell(nJobsW,1);
    for i = 1:nJobsW
        k = work2Job(i);
        cik = callArroba(ciW,{k});
        
        objW{i} = arroba(fW,...
                         [1,2,3],...
                         {...
                         objJk,...
                         stepSchedulesW(i),...
                         wellSol,...
                        'xScale',...
                         xScale,...
                         'vScale',...
                         vScale,...
                         'uScale',...
                         cellControlScales{cik}...
                         },true);
    end
    obj = objW;
end
targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,jobSchedule,work2Job,varargin{:});


%%% objective function on the client side (for plotting!)
objClient = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    objClient{k} = arroba(@mrstTimePointFuncWrapper,...
        [1,2,3],...
        {...
        objJk,...
        stepSchedules(k),...
        wellSol,...
        'xScale',...
        xScale,...
        'vScale',...
        vScale,...
        'uScale',...
        cellControlScales{callArroba(ci,{k})}...
        },true);
end


%%  Bounds for all variables!
maxProd = struct('ORAT',500*meter^3/day,'WRAT',500*meter^3/day,'LRAT',500*meter^3/day,'RESV',inf,'BHP',420*barsa);
minProd = struct('ORAT',eps,  'WRAT',eps,  'LRAT',eps,  'RESV',eps,  'BHP',(380)*barsa);
maxInj = struct('RATE',500*meter^3/day,'RESV',inf,  'BHP',410*barsa);
minInj = struct('RATE',eps,'RESV',eps,  'BHP',380*barsa);


% Control input bounds for all wells!
[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProd,'minProd',minProd,...
    'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);
lbu = schedules2CellControls(lbSchedules,'doScale',true,'cellControlScales',cellControlScales);
ubu = schedules2CellControls(ubSchedules,'doScale',true,'cellControlScales',cellControlScales);

% wellSol bounds  (Algebraic variables bounds)
[ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,maxProd,maxInj,minProd,minInj);
ubvS = wellSol2algVar(ubWellSol,'doScale',true,'vScale',vScale);
lbvS = wellSol2algVar(lbWellSol,'doScale',true,'vScale',vScale);
lbv = repmat({lbvS},totalPredictionSteps,1);
ubv = repmat({ubvS},totalPredictionSteps,1);

%%%%%%%%%%% Initialization lower and upper - bounds
maxState = struct('pressure',(400)*barsa,'s',1);
minState = struct('pressure',(380)*barsa,'s',0.0999);                     
[ubState] = stateBounds(reservoirP.state,maxState);
[lbState] = stateBounds(reservoirP.state,minState);
lbxS = stateMrst2stateVector( lbState,'doScale',true,'xScale',xScale );
ubxS = stateMrst2stateVector( ubState,'doScale',true,'xScale',xScale );
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Max saturation for well producer cells

% get well perforation grid block set
prodInx  = (vertcat(wellSol.sign) < 0);
wc    = vertcat(W(prodInx).cells);

maxSat = struct('pressure',inf,'s',0.65); 
[satWMax] = stateBounds(ubState,maxSat,'cells',wc);
ubxS = stateMrst2stateVector( satWMax,'doScale',true,'xScale',xScale );
ubxsatWMax = repmat({ubxS},totalPredictionSteps,1);
ubx = cellfun(@(x1,x2)min(x1,x2),ubxsatWMax,ubx,'UniformOutput',false);



%% Initial Active set!
initializeActiveSet = true;
if initializeActiveSet
[ lowActive,upActive ] = activeSetFromWells( reservoirP,totalPredictionSteps);
else
    lowActive = [];
    upActive = [];
end





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
simFunc =@(sch) runScheduleADI(reservoirP.state, reservoirP.G, reservoirP.rock, reservoirP.system, sch);


wc    = vertcat(W.cells);
fPlot = @(x)[max(x);min(x);x(wc)];

%prodInx  = (vertcat(wellSol.sign) < 0);
%wc    = vertcat(W(prodInx).cells);
%fPlot = @(x)x(wc);

plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d,ss,objClient,times,xScale,cellControlScales,vScale,cellControlScalesPlot,controlSchedules,wellSol,ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true,'plotSchedules',false,'pF',fPlot,'sF',fPlot,varargin{:});

%%  Initialize from previous solution?

if exist('optimalVars.mat','file') == 2
    load('optimalVars.mat','x','u','v');
elseif exist('itVars.mat','file') == 2
    load('itVars.mat','x','u','v');
else
	x = [];
    v = [];
    u  = schedules2CellControls( controlSchedules,'doScale',true,'cellControlScales',cellControlScales);
    %[x] = repmat({ss.state},totalPredictionSteps,1);
end


%% call REMSO

[u,x,v,f,d,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
    'tol',1e-2,'lkMax',4,'debugLS',true,...
    'lowActive',lowActive,'upActive',upActive,...
    'plotFunc',plotSol,'max_iter',500,'x',x,'v',v,'debugLS',false,'saveIt',true);

%% plotSolution
plotSol(x,u,v,d,'simFlag',true);

%% Save result?
%save optimalVars u x v f d M
