% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global

% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi ad-core

% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/procedural'));
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/robust'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));

%% Initialize reservoir -  the Simple reservoir
[reservoirP] = initReservoir( 'simple10x1x10.data','Verbose',true);
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


%%% The sum of the last elements in the algebraic variables is the objective
nCells = reservoirP.G.cells.num;
stepNPV = arroba(@NPVStepM,[1,2],{nCells,'scale',1/10000,'sign',-1},true);
vScale = [vScale;1];

%% Instantiate the simulators for each interval, locally and for each worker.

% ss.stepClient = local (client) simulator instances
% ss.state = scaled initial state
% ss.nv = number of algebraic variables
% ss.ci = piecewice control mapping on the client side
% ss.step =  worker simulator instances

nR = 10;
ss = cell(nR,1);
for r = 1:nR
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
    
    
    ss{r}.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
    
    ss{r}.state = ss{r}.state + rand(size(ss{r}.state))/100;  %% add uncertainty
    
    ss{r}.step = step;
    ss{r}.ci = ci;
    ss{r}.outputF = @npvStages;
end

obj = @objSumRisks;

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
maxProd = struct('ORAT',inf,'WRAT',inf,'GRAT',inf,'BHP',inf);
minProd = struct('ORAT',-inf,  'WRAT',-inf,  'GRAT',-inf,'BHP',-inf);
maxInj = struct('ORAT',inf,'WRAT',inf,'GRAT',inf,'BHP',inf);
minInj = struct('ORAT',-inf,  'WRAT',-inf,  'GRAT',-inf,'BHP',-inf);

% wellSol bounds  (Algebraic variables bounds)
[ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
    'maxProd',maxProd,...
    'maxInj',maxInj,...
    'minProd',minProd,...
    'minInj',minInj);
ubvS = wellSol2algVar(ubWellSol,'vScale',vScale);
lbvS = wellSol2algVar(lbWellSol,'vScale',vScale);
lbv = repmat({[lbvS;-inf]},totalPredictionSteps,1);
ubv = repmat({[ubvS;-inf]},totalPredictionSteps,1);

% State lower and upper - bounds
maxState = struct('pressure',inf,'sW',inf);
minState = struct('pressure',-inf,'sW',-inf);
ubxS = setStateValues(maxState,'nCells',nCells,'xScale',xScale);
lbxS = setStateValues(minState,'nCells',nCells,'xScale',xScale);
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);


%%  Initialize from previous solution?

if exist('optimalVars.mat','file') == 2
    load('optimalVars.mat','x','u','v');
elseif exist('itVars.mat','file') == 2
    load('itVars.mat','x','u','v');
else
    u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);
end


objectiveSS = @(u,varargin) simulateSystemSS_R(u,ss,obj,varargin{:});


uDim = cellfun(@(x)size(x,1),u);
[outputCons,lbC,ubC,consSparsity] = outputVarsBoundSelector(lbx,ubx,lbv,ubv,uDim,ci);

consSizes = cellfun(@(x)size(x,1),lbC);
cons = cell(numel(consSizes),1);
for k=1:size(cons)
    cons{k} = @(xsk,vsk,uk,varargin) concatenateTargetK(k,xsk,vsk,uk,outputCons{k},consSizes,varargin{:});
end

constraintSS = @(u,varargin) simulateSystemSS_R(u,ss,cons,varargin{:});


x0         = cell2mat(u);   % The starting point.
options.lb = cell2mat(lbu);  % Lower bound on the variables.
options.ub = cell2mat(ubu);  % Upper bound on the variables.
%options.cl = cell2mat(lbC);   % Lower bounds on the constraint functions.
%options.cu = cell2mat(ubC);   % Upper bounds on the constraint functions.

[ fM ] = memorizeLastSimulation(u,[],true);


fwdObj = @(x) ssFwdMemory(mat2cell(x,uDim,1),...
    @(xx,varargin)objectiveSS(xx,'gradients',false,varargin{:}),...
    fM,...
    'replace',false);
gradObj = @(x) cell2mat(ssFwdMemory(mat2cell(x,uDim,1),...
    @(xx,varargin)objectiveSS(xx,'gradients',true,varargin{:}),...
    fM,...
    'replace',true,'gradFirst',true));
fwdCons = @(x) ssFwdMemory(mat2cell(x,uDim,1),...
    @(xx,varargin)constraintSS(xx,'gradients',false,varargin{:}),...
    fM,...
    'replace',false);
gradCons = @(x) sparse(cell2mat(ssFwdMemory(mat2cell(x,uDim,1),...
    @(xx,varargin)constraintSS(xx,'gradients',true,varargin{:}),...
    fM,...
    'replace',true,'gradFirst',true)));

% The callback functions.
funcs.objective        = fwdObj;
funcs.gradient         = gradObj;
%funcs.constraints       = fwdCons;
%funcs.jacobian          = gradCons;
%funcs.jacobianstructure = @(x) sparse(cell2mat(consSparsity));
%funcs.iterfunc         = @callback;

% Set the IPOPT options.
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.tol         = 1e-7;
options.ipopt.max_iter    = 100;

% Run IPOPT.
[x info] = ipopt(x0,funcs,options);

