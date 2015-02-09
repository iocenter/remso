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
initPool('restart',true);


%% Initialize reservoir the Egg model
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
ci  = arroba(@controlIncidence,2,{reservoirP.schedule.step.control});


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


%%% The sum of the last elements in the algebraic variables is the objective
nCells = reservoirP.G.cells.num;
stepNPV = arroba(@NPVStepM,[1,2],{nCells,'scale',1/100000,'sign',-1},true);

vScale = [vScale;1];

%% Instantiate the simulators for each interval, locally and for each worker.

% ss.stepClient = local (client) simulator instances 
% ss.state = scaled initial state
% ss.nv = number of algebraic variables
% ss.ci = piecewice control mapping on the client side
% ss.jobSchedule = Step distribution among workers client side
% ss.work2Job = Step distribution among workers worker side
% ss.step =  worker simulator instances 

stepClient = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    ss.stepClient{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
                                        'xScale',xScale,...
                                        'vScale',vScale,...
                                        'uScale',cellControlScales{cik},...
                                        'algFun',stepNPV,...
										'saveJacobians',false,...
                                        varargin{:});
end



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
		    'algFun',stepNPV,...
            'saveJacobians',false...
            },...
            true);
    end
    step = stepW;
end

ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
ss.jobSchedule = jobSchedule;
ss.work2Job = work2Job;
ss.step = step;
ss.ci = ci;




%% instantiate the objective function


%%% objective function on the client side (for plotting!)
objClient = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    objClient{k} = arroba(@lastAlg,[1,2,3],{},true);
end
%%% Investigate if it is efficient to evaluate the objective in the workers
spmd
    nJobsW = numel(work2Job);
    objW = cell(nJobsW,1);
    for i = 1:nJobsW
        objW{i} = arroba(@lastAlg,[1,2,3],{},true);
    end
    obj = objW;
end
targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,jobSchedule,work2Job,varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hard constraints REMSO only


%% Bounds for all wells!
maxProdH = struct('ORAT',inf*meter^3/day,'WRAT',inf*meter^3/day,'GRAT',inf*meter^3/day,'BHP',inf*barsa);
minProdH = struct('ORAT',-inf*meter^3/day,  'WRAT',-inf*meter^3/day,  'GRAT',-inf*meter^3/day,'BHP',-inf*barsa);
maxInjH = struct('ORAT',inf*meter^3/day,'WRAT',inf*meter^3/day,'GRAT',inf*meter^3/day,'BHP',inf*barsa);
minInjH = struct('ORAT',-inf*meter^3/day,  'WRAT',-inf*meter^3/day,  'GRAT',-inf*meter^3/day,'BHP',-inf*barsa);

% wellSol bounds  (Algebraic variables bounds)
[ubWellSolH,lbWellSolH] = wellSolScheduleBounds(wellSol,...
    'maxProd',maxProdH,...
    'maxInj',maxInjH,...
    'minProd',minProdH,...
    'minInj',minInjH);
ubvSH = wellSol2algVar(ubWellSolH,'vScale',vScale);
lbvSH = wellSol2algVar(lbWellSolH,'vScale',vScale);
lbvH = repmat({lbvSH},totalPredictionSteps,1);
ubvH = repmat({ubvSH},totalPredictionSteps,1);

% State lower and upper - bounds
maxStateH = struct('pressure',inf*barsa,'sW',1);
minStateH = struct('pressure',0*barsa,'sW',0);
ubxSH = setStateValues(maxStateH,'nCells',nCells,'xScale',xScale);
lbxSH = setStateValues(minStateH,'nCells',nCells,'xScale',xScale);
lbxH = repmat({[lbxSH;-inf]},totalPredictionSteps,1);
ubxH = repmat({[ubxSH;inf]},totalPredictionSteps,1);


%%  Bounds for all variables!
maxProdInput = struct('BHP',420*barsa);
minProdInput = struct('BHP',(380)*barsa);
maxInjInput = struct('RATE',500*meter^3/day);
minInjInput = struct('RATE',2*meter^3/day);


% Control input bounds for all wells!
[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProdInput,'minProd',minProdInput,...
    'maxInj',maxInjInput,'minInj',minInjInput,'useScheduleLims',false);
lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales);

% Bounds for all wells!
maxProd = struct('ORAT',500*meter^3/day,'WRAT',500*meter^3/day,'GRAT',500*meter^3/day,'BHP',420*barsa);
minProd = struct('ORAT',2*meter^3/day,  'WRAT',0*meter^3/day,  'GRAT',-inf*meter^3/day,'BHP',380*barsa);
maxInj = struct('ORAT',500*meter^3/day,'WRAT',500*meter^3/day,'GRAT',500*meter^3/day,'BHP',410*barsa);
minInj = struct('ORAT',0*meter^3/day,  'WRAT',2*meter^3/day,  'GRAT',-inf*meter^3/day,'BHP',380*barsa);

% wellSol bounds  (Algebraic variables bounds)
[ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
    'maxProd',maxProd,...
    'maxInj',maxInj,...
    'minProd',minProd,...
    'minInj',minInj);
ubvS = wellSol2algVar(ubWellSol,'vScale',vScale);
lbvS = wellSol2algVar(lbWellSol,'vScale',vScale);
lbv = repmat({[lbvS;-inf]},totalPredictionSteps,1);
ubv = repmat({[ubvS;inf]},totalPredictionSteps,1);

%%%%%%%%%%% Initialization lower and upper - bounds
maxState = struct('pressure',400*barsa,'sW',1);
minState = struct('pressure',380*barsa,'sW',0.0999);
lbxS = setStateValues( minState,'nCells',nCells,'xScale',xScale);
ubxS = setStateValues( maxState,'nCells',nCells,'xScale',xScale);
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Max saturation for well producer cells

% get well perforation grid block set
prodInx  = (vertcat(wellSol.sign) < 0);
wc    = vertcat(W(prodInx).cells);

maxSat = struct('pressure',inf,'sW',0.65); 
ubxS = setStateValues(maxSat,'x',ubxS,'xScale',xScale,'cells',wc,'nCells',reservoirP.G.cells.num);
ubxsatWMax = repmat({ubxS},totalPredictionSteps,1);
ubx = cellfun(@(x1,x2)min(x1,x2),ubxsatWMax,ubx,'UniformOutput',false);



%% Initial Active set!
initializeActiveSet = true;
if initializeActiveSet
    vDims = cellfun(@numel,lbv);
    [ lowActive,upActive ] = activeSetFromWells(vDims,reservoirP,totalPredictionSteps);
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
simFunc =@(sch,varargin) runScheduleADI(reservoirP.state, reservoirP.G, reservoirP.rock, reservoirP.system, sch,'force_step',false,varargin{:});


wc    = vertcat(W.cells);
fPlot = @(x)[max(x);min(x);x(wc)];

%prodInx  = (vertcat(wellSol.sign) < 0);
%wc    = vertcat(W(prodInx).cells);
%fPlot = @(x)x(wc);

plotSol = @(x,u,v,xd,varargin) plotSolution( x,u,v,xd,ss,objClient,times,xScale,cellControlScales,vScale,cellControlScalesPlot,controlSchedules,wellSol,ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true,'plotSchedules',false,'pF',fPlot,'sF',fPlot,varargin{:});

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

%  Exploit a bit more of structure, include input bounds to the reservoir
%  states too!

maxStateI = struct('pressure',inf*barsa,'sW',0.95);
minStateI = struct('pressure',minProdInput.BHP,'sW',0.05);
lbxSI = setStateValues( minStateI,'nCells',nCells,'xScale',xScale);
ubxSI = setStateValues( maxStateI,'nCells',nCells,'xScale',xScale);
lbxI = repmat({lbxSI},totalPredictionSteps,1);
ubxI = repmat({ubxSI},totalPredictionSteps,1);


x0 = stateMrst2stateVector( reservoirP.state,'xScale',xScale);  % initial state must be feasible!
lbx = cellfun(@(x1,x2)min(max(x1,x2),x0),lbxI,lbx,'UniformOutput',false);
ubx = cellfun(@(x1,x2)max(min(x1,x2),x0),ubxI,ubx,'UniformOutput',false);
[u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
                                             'lbxH',lbxH,'ubxH',ubxH,'lbvH',lbvH,'ubvH',ubvH,...
    'tol',1e-2,'lkMax',4,'debugLS',true,...
    'lowActive',lowActive,'upActive',upActive,...
    'plotFunc',plotSol,'max_iter',500,'x',x,'v',v,'debugLS',false,'saveIt',true,'condensingParallel',false);

%% plotSolution
plotSol(x,u,v,xd,'simFlag',true);

%% Save result?
%save optimalVars u x v f xd M
