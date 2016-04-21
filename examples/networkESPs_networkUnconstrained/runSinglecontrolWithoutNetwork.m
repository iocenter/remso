% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
clc
clear
clear global

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
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/remsoCrossSequential'));
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));


%% Initialize reservoir -  the Simple reservoir
[reservoirP] = initReservoir('RATE10x5x10.txt', 'Verbose',true);

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
ci  = arroba(@controlIncidence, 2 ,{reservoirP.schedule.step.control});

%% Variables Scaling
xScale = setStateValues(struct('pressure',5*barsa,'sW',0.01),'nCells',nCells);

if (isfield(reservoirP.schedule.control,'W'))
    W =  reservoirP.schedule.control.W;
else
    W = processWells(reservoirP.G, reservoirP.rock,reservoirP.schedule.control(1),'DepthReorder', true);
end
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

freqScale = [15;15;15;15;15]; % in Hz
flowScale = [5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day);5*(meter^3/day); ...
    5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day);5*(meter^3/day)];

pressureScale = [5*barsa;5*barsa;5*barsa;5*barsa;5*barsa];

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
stepNPV = arroba(@NPVStepM,[1,2, 3],{nCells,'scale',1/100000,'sign',-1},true);

vScale = [vScale; 1];

step = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
        'xScale',xScale,...
        'vScale',vScale,...
        'uScale',cellControlScales{cik},...
        'algFun',stepNPV,...
        'fixedWells', fixedWells, ...
        'saveTargetJac', true,...
        varargin{:});
end

ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
ss.step = step;
ss.ci = ci;

%% instantiate the objective function
obj = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    obj{k} = arroba(@lastAlg,[1,2,3],{},true);
end
targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,varargin{:});

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
lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);


% Bounds for all wells!
minProd = struct('ORAT', 5*meter^3/day,  'WRAT', 0*meter^3/day,  'GRAT', -inf,'BHP',100*barsa);
maxProd = struct('ORAT',200*meter^3/day,'WRAT',200*meter^3/day,'GRAT', inf,'BHP',500*barsa);

minInj = struct('ORAT',-inf,  'WRAT', 5*meter^3/day,  'GRAT', -inf,'BHP', 5*barsa);
maxInj = struct('ORAT',inf,'WRAT', 300*meter^3/day,'GRAT', inf,'BHP',800*barsa);

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

plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d, lbv, ubv, lbu, ubu, ss,obj,times,xScale,cellControlScales,vScale, nScale, ...
    cellControlScalesPlot,controlSchedules,wellSol, netSol, ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true, 'plotNetsol', true, ...
    'numNetConstraints', numel(nScale), 'plotNetControls', false, 'numNetControls', numel(pScale), 'freqCst', numel(freqScale), 'pressureCst',numel(pressureScale),  'flowCst',numel(flowScale), ...
    'plotSchedules',false,'pF',fPlot,'sF',fPlot, 'fixedWells', fixedWells, 'extremePoints', extremePoints, 'plotCumulativeObjective', true, 'qlMin', qlMin,  'qlMax', qlMax, 'nStages', numStages, ...
    'freqMin', freqMin, 'freqMax', freqMax, 'baseFreq', baseFreq, 'reservoirP', reservoirP, 'plotNetwork', true, 'dpFunction', @dpBeggsBrillJDJ,  varargin{:});

% remove network control to initialize well controls vector (w)
cellControlScales = cellfun(@(w) w(1:end-numel(p)) ,cellControlScales, 'UniformOutput', false);

%%  Initialize from previous solution?
x = [];
v = [];
u = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);

testFlag = false;
if testFlag
    addpath(genpath('../../optimization/testFunctions'));
    [~, ~, ~, simVars, xs, vs] = simulateSystemSS(u, ss, []);
    [ei, fi, vi] = testProfileGradients(xs,u,vs,ss.step,ss.ci,ss.state, 'd', 1, 'pert', 1e-5, 'all', false);
end

optmize = true;
loadPrevSolution = false;
plotSolution = false;

algorithm = 'remso';
switch algorithm
    case 'remso'
        if loadPrevSolution
            load itVars;
        end
        
        if optmize
            %% call REMSO
            [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
                'tol',1e-6,'lkMax',4,'debugLS',true,...
                'lowActive',lowActive,'upActive',upActive,...
                'plotFunc',plotSol,'max_iter', 500,'x',x,'v',v,'debugLS',false,'saveIt',true, 'computeCrossTerm', false, 'condense', true);
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
        
    case 'greedyGeneral'
        if isempty(x)
            x = cell(totalPredictionSteps,1);
        end
        if isempty(v)
            v = cell(totalPredictionSteps,1);
        end
        xd = cell(totalPredictionSteps,1);
        ssK = ss;
        
        uK = u(1);
        kFirst = 1;
        iC = 1;
        if loadPrevSolution
            load greedyStrategy;
        end
        
        recoverPreviousSolution = false;
        if recoverPreviousSolution
            load greedyStrategy.mat;
            iC = kLast;
            lastControlSteps = lastControlSteps(kLast:end);
            kFirst = kLast;
            ssK.state = lastState;
        end
        
        if optimize
            for kLast = lastControlSteps'
                kLast
                if loadPrevSolution
                    uK = u(iC);
                    xK = x(kFirst:kLast);
                    vK = v(kFirst:kLast);
                else
                    xK = [];
                    vK = [];
                    if iC > 1
                        uK = u(iC-1);
                    end
                end
                
                totalPredictionStepsK = kLast-kFirst+1;
                
                ssK.step = ss.step(kFirst:kLast);
                ssK.ci  = arroba(@controlIncidence,2,{ones(totalPredictionStepsK,1)});
                
                lbxK = lbx(kFirst:kLast);
                ubxK = ubx(kFirst:kLast);
                lbvK = lbv(kFirst:kLast);
                ubvK = ubv(kFirst:kLast);
                lbuK = lbu(iC);
                ubuK = ubu(iC);
                
                obj = cell(totalPredictionStepsK,1);
                for k = 1:totalPredictionStepsK
                    obj{k} = arroba(@lastAlg,[1,2,3],{},true);
                end
                targetObjK = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ssK,varargin{:});
                
                
                [ukS,xkS,vkS,f,xdK,M,simVars,converged] = remso(uK,ssK,targetObjK,'lbx',lbxK,'ubx',ubxK,'lbv',lbvK,'ubv',ubvK,'lbu',lbuK,'ubu',ubuK,...
                    'tol',1e-6,'lkMax',4,'debugLS',false,...
                    'skipRelaxRatio',inf,...
                    'lowActive',[],'upActive',[],...
                    'plotFunc',plotSol,'max_iter', 500,'x',xK,'v',vK,'saveIt',false, 'condense', true,'computeCrossTerm',false, 'qpAlgorithm', 1);
                
                x(kFirst:kLast) = xkS;
                v(kFirst:kLast) = vkS;
                xd(kFirst:kLast) = xdK;
                u(iC) = ukS;
                
                kFirst = kLast+1;
                ssK.state = xkS{end};
                lastState = ssK.state;
                iC = iC+1;
                save('greedyStrategy.mat', 'lastState', 'u', 'kLast');
            end
            save('greedyStrategy.mat','x', 'xd', 'v', 'u');
        end
        
        if plotSolution
            if ~optimize
                [~, ~, ~, simVars, x, v] = simulateSystemSS(u, ss, []);
            end
            xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
            plotSol(x,u,v,xd, 'simFlag', false);
        end
    otherwise
        error('algorithm must be either remso or greedyGeneral')
        
end