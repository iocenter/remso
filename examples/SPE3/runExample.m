% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.
 
% Make sure the workspace is clean before we start
clc
clear
clear global

% Required MRST modules
mrstModule clear
mrstModule add deckformat
mrstModule add ad-fi ad-core ad-blackoil ad-props

% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/procedural'));
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/remsoCrossSequential'));
addpath(genpath('../../optimization/testFunctions'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));


[reservoirP] = initReservoir();


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

%
nCells = reservoirP.G.cells.num;

%% Variables Scaling

xScale = setStateValues(struct('pressure',100*psia,...
                               'sW',0.01,...
                               'rGH',0.01),...
                        'nCells',nCells);


W =  reservoirP.schedule.control.W;
wellSol = initWellSolLocal(W, reservoirP.state);
vScale = wellSol2algVar( wellSolScaling(wellSol,'activeComponents',reservoirP.system.activeComponents,...
    'bhp',100*psia,...
    'qWs',10*stb/day,...
    'qOs',100*stb/day,...
    'qGs',100*(10*ft)^3/day),...
    'activeComponents',reservoirP.system.activeComponents);

cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
    'RATE',100*(10*ft)^3/day,...  % injector contrlled by gas rate
    'GRAT',100*(10*ft)^3/day,...    
    'ORAT',100*stb/day,...
    'WRAT',10*stb/day,...
    'LRAT',0,...
    'RESV',0,...
    'BHP',100*psia));

%% instantiate the objective function as an aditional Algebraic variable
%%% The sum of the last elements in the algebraic variables is the objective
nCells = reservoirP.G.cells.num;
stepNPV = arroba(@NPVVOStep,[1,2],{nCells,'scale',1e-8,'sign',-1,'OilPrice',300},true);
%% Instantiate the simulators for each interval, locally and for each worker.

% ss.stepClient = local (client) simulator instances
% ss.state = scaled initial state
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



ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale,...
    'activeComponents',reservoirP.system.activeComponents,...
    'fluid',reservoirP.fluid);
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

% Bounds for Controls!   
maxProd = struct('GRAT',6200*(1000 * ft^3)/day);
minProd = struct('GRAT',0*(1000 * ft^3)/day);
maxInj = struct('RATE',4700*(1000 * ft^3)/day); % Mscf
minInj = struct('RATE',0*(1000 * ft^3)/day);


% Control input bounds for all wells!
[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProd,'minProd',minProd,...
    'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);
lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales);

% Bounds for all wells!
maxProd = struct('ORAT',inf,'WRAT',inf,'GRAT',inf,'BHP',inf);
minProd = struct('ORAT',0,  'WRAT',0,  'GRAT',0,'BHP',1050*psia);
maxInj = struct('ORAT',inf,'WRAT',inf,'GRAT',inf,'BHP',4000*psia); 
minInj = struct('ORAT',0,  'WRAT',0,  'GRAT',0,'BHP',0);

ubv = cell(totalPredictionSteps,1);
lbv = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    % during a shooting period the number of wells cannot change, sorry...
    W = stepSchedules(k).control(1).W;
    wellSol = initWellSolLocal(W, reservoirP.state);

% wellSol bounds  (Algebraic variables bounds)
    [ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
        'maxProd',maxProd,...
        'maxInj',maxInj,...
        'minProd',minProd,...
        'minInj',minInj);
    
    ubwk = wellSol2algVar(ubWellSol,'vScale',vScale,'activeComponents',reservoirP.system.activeComponents);
    lbwk = wellSol2algVar(lbWellSol,'vScale',vScale,'activeComponents',reservoirP.system.activeComponents);

    %%%%  
    % Implement here the bound for the bounds for algFun at each shooting period
    % in this example we only have the objective (which is unbounded!)
    
    ubn = inf;
    lbn = -inf;
    
    lbv{k} = [lbwk;lbn];
    ubv{k} = [ubwk;ubn];
    
    
end
% State lower and upper - bounds
maxState = struct('pressure',3600*psia,...
                  'sW',1-2*sqrt(eps),...
                  'rGH',1);
minState = struct('pressure',1215*psia,...
                  'sW',0,...
                  'rGH',0);
lbxS = setStateValues(minState,...
                      'nCells',nCells,...
                      'xScale',xScale);
ubxS = setStateValues(maxState,...
                      'nCells',nCells,...
                      'xScale',xScale);
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);




%% Initial Active set!
initializeActiveSet = false;
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

plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d,ss,obj,times,xScale,cellControlScales,vScale,cellControlScalesPlot,controlSchedules,wellSol,ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true,'plotSchedules',false,'pF',fPlot,'sF',fPlot,...
    'activeComponents',reservoirP.system.activeComponents,...
    'fluid',reservoirP.fluid,...
    varargin{:});

%%  Initialize from previous solution?

if exist('optimalVars.mat','file') == 2
    load('optimalVars.mat','x','u','v','xd');
elseif exist('itVars.mat','file') == 2
    load('itVars.mat','x','u','v','xd');
else
	x = [];
    v = [];
    u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);
    %[x] = repmat({ss.state},totalPredictionSteps,1);
end


%% call REMSO

[u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbxH',lbx,'ubxH',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
    'tol',1e-6,'lkMax',4,'debugLS',true,'condT',1e9,...
    'lowActive',lowActive,'upActive',upActive,...
    'plotFunc',plotSol,'plot',false,'max_iter',500,'x',x,'v',v,'debugLS',false,'saveIt',true,'condense',true,'computeCrossTerm',false,'nCons',10);

%% plotSolution
plotSol(x,u,v,xd,'simFlag',true);

%% Save result?
%save optimalVars u x v f xd M
