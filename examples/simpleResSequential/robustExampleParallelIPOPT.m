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
addpath(genpath('../../optimization/parallel'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/robust'));
addpath(genpath('../../optimization/robustParallel'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));

mrstVerbose off;

% Open a matlab pool depending on the machine availability
initPool('restart',true);

nR = 8;
%%  Who will do what - Distribute the computational effort!
nWorkers = getNumWorkers();
if nWorkers == 0
    nWorkers = 1;
end
[ jobSchedule ] = divideJobsSequentially(nR ,nWorkers);
jobSchedule.nW = nWorkers;

work2Job = distributeVariables(jobSchedule.work2Job,jobSchedule);

pScale = 5*barsa;
sScale = 0.01;
qwScale = 10*meter^3/day;
qoScale = 10*meter^3/day;
objScale = 1/100000;

% to get initial schedule only
[reservoirP] = initReservoir( 'simple10x1x10.data','Verbose',true);
schedule = reservoirP.schedule;
clear reservoirP

lastControlSteps = findControlFinalSteps( schedule.step.control );
controlSchedules = multipleSchedules(schedule,lastControlSteps);

cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
    'RATE',qwScale,...
    'ORAT',qoScale,...
    'WRAT',qwScale,...
    'LRAT',qwScale,...
    'RESV',0,...
    'BHP',pScale));

spmd

	ss = cell(numel(work2Job),1);

    for r=1:numel(work2Job)
    
    % read the realization
    [reservoirP] = initReservoir( 'simple10x1x10.data','Verbose',true);
    nCells = reservoirP.G.cells.num;
    
    %% Multiple shooting problem set up
    totalPredictionSteps = numel(reservoirP.schedule.step.val);  % MS intervals
    stepSchedules = multipleSchedules(reservoirP.schedule,1:totalPredictionSteps);
    
    
    % Piecewise linear control -- mapping the step index to the corresponding
    % control
    ci  = arroba(@controlIncidence,2,{reservoirP.schedule.step.control});
    
    
    %% Variables Scaling
    xScale = setStateValues(struct('pressure',pScale,'sW',sScale),'nCells',nCells);
    
    
    W =  reservoirP.schedule.control.W;
    wellSol = initWellSolLocal(W, reservoirP.state);
    wellSolScales = wellSol2algVar( wellSolScaling(wellSol,'bhp',pScale,'qWs',qwScale,'qOs',qoScale) );
    
    stepNPV = arroba(@NPVStepM,[1,2],{nCells,'scale',objScale,'sign',-1},true);

    
    nW = numel(W);
    vflow = arroba(@averageFlowRate,[1,2],{nCells,'scale',-1/(qwScale)},true);
    
    algFun = concatenateMrstTargets([vflow,stepNPV],false,[nW,1]);
    
    extraAlgScales = ones(nW+1,1);
    
    vScale = [wellSolScales;extraAlgScales];
    
    
    step = cell(totalPredictionSteps,1);
    for k=1:totalPredictionSteps
        
        cik = callArroba(ci,{k});
        step{k} = arroba(@mrstStep,...
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
            'algFun',algFun,...
            'saveJacobians',false...
            },...
            true);
    end
    
    
    sW = 0.2 + 0.3 * rand(nCells,1);
    reservoirP.state.s = [sW,1-sW];
    
    ss{r}.step = step;
    ss{r}.ci = ci;
    ss{r}.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
    ss{r}.outputF = arroba(@lastNV,[1,2,3],{nW+1},true);
    end  % r=1:numel(work2Job)

end % spmd
sss.ss = ss;
sss.nR = nR;
sss.jobSchedule = jobSchedule;
sss.eta = 0.8;

nW = numel(schedule.control(1).W);
totalPredictionSteps = numel(schedule.step.val);
selection = repmat([false(nW,1);true],totalPredictionSteps,1);
obj = @(s,u,varargin)sumSelectionS(s,u,selection,varargin{:});



%%  Bounds for all variables!
maxProdInput = struct('BHP',200*barsa);
minProdInput = struct('BHP',(50)*barsa);
maxInjInput = struct('RATE',250*meter^3/day);
minInjInput = struct('RATE',0);


% Control input bounds for all wells!
[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProdInput,'minProd',minProdInput,...
    'maxInj',maxInjInput,'minInj',minInjInput,'useScheduleLims',false);
lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales);

totalPredictionSteps = numel(schedule.step.val);
lbs = repmat(-inf(nW+1,1),totalPredictionSteps,1);
ubs = repmat([zeros(nW,1);inf],totalPredictionSteps,1);


u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);


objectiveSS = @(u,varargin) simulateSystemSS_R(u,sss,obj,varargin{:});


uDims = cellfun(@numel,u);
[outputCons,lbC,ubC,consSparsity] = stochasticOutputVarsSelector(lbs,ubs,uDims);



constraintSS = @(u,varargin) simulateSystemSS_R(u,sss,outputCons,varargin{:});


x0         = cell2mat(u);   % The starting point.
options.lb = cell2mat(lbu);  % Lower bound on the variables.
options.ub = cell2mat(ubu);  % Upper bound on the variables.
options.cl = lbC;   % Lower bounds on the constraint functions.
options.cu = ubC;   % Upper bounds on the constraint functions.

[ fM ] = memorizeLastSimulation(u,[],true);


fwdObj = @(x) ssFwdMemory(mat2cell(x,uDims,1),...
    @(xx,varargin)objectiveSS(xx,'gradients',false,varargin{:}),...
    fM,...
    'replace',false);
gradObj = @(x) cell2mat(ssFwdMemory(mat2cell(x,uDims,1),...
    @(xx,varargin)objectiveSS(xx,'gradients',true,varargin{:}),...
    fM,...
    'replace',true,'gradFirst',true));
fwdCons = @(x) ssFwdMemory(mat2cell(x,uDims,1),...
    @(xx,varargin)constraintSS(xx,'gradients',false,varargin{:}),...
    fM,...
    'replace',false);
gradCons = @(x) sparse(cell2mat(ssFwdMemory(mat2cell(x,uDims,1),...
    @(xx,varargin)constraintSS(xx,'gradients',true,varargin{:}),...
    fM,...
    'replace',true,'gradFirst',true)));

% The callback functions.
funcs.objective        = fwdObj;
funcs.gradient         = gradObj;
funcs.constraints       = fwdCons;
funcs.jacobian          = gradCons;
funcs.jacobianstructure = @(x) sparse(consSparsity);
%funcs.iterfunc         = @callback;

% Set the IPOPT options.
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.tol         = 1e-7;
options.ipopt.max_iter    = 100;

% Run IPOPT.
[x info] = ipopt(x0,funcs,options);

