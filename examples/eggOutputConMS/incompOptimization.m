% Make sure the workspace is clean before we start
clc
clear
clear global

% Required MRST modules
mrstModule add deckformat
mrstModule add mimetic adjoint

% Include REMSO functionalities
% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstDerivated/modules/adjoint'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/adjoint'));
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/remsoCrossSequential'));
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));



mrstVerbose off;

% Open a matlab pool depending on the machine availability


[reservoirP] = loadIncompressibleEgg('./reservoirData/');


% do not display reservoir simulation information!
mrstVerbose off;


%% Multiple shooting problem set up
totalPredictionSteps = numel(reservoirP.schedule);  % MS intervals



stepSchedules = reservoirP.schedule;


% Piecewise linear control -- mapping the step index to the corresponding
% control
ci  = arroba(@controlIncidence,2,{reservoirP.simulationControlSteps});


qScale = 10*meter^3/day;
pScale = 5*barsa;
sScale = 0.01;


xScale = sScale*ones(reservoirP.G.cells.num,1);
vScale = qScale*ones(numel(vertcat(reservoirP.W.cells)),1);
uScale = ones(numel(reservoirP.W),1);


for ui = 1:numel(uScale)
    if strcmp(reservoirP.W(ui).type,'bhp')
        uScale(ui) = pScale;
    elseif strcmp(reservoirP.W(ui).type,'rate')
        uScale(ui) = qScale;
    else
        error('what kind of control?')
    end
end



schedule = reservoirP.schedule;
controls = initControls(schedule);

stepNPV = arroba(@simpleNPVADI,[6,7],{reservoirP.G, reservoirP.S, reservoirP.W, reservoirP.rock, reservoirP.fluid, controls,'scale',1e-7,'sign',-1},true);
vScale = [vScale;1];



step = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,reservoirP.state.wellSol,stepSchedules(k),reservoirP,...
                                        'xScale',xScale,...
                                        'vScale',vScale,...
                                        'uScale',uScale,...
                                        'algFun',stepNPV,...
                                        varargin{:});
end



ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
ss.step = step;
ss.ci = ci;


%%% objective function
obj = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    obj{k} = arroba(@lastAlg,[1,2,3],{},true);
end
targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,varargin{:});


injIndex = vertcat(reservoirP.W.sign) == 1;
ProdIndex = vertcat(reservoirP.W.sign) == -1;


nInj = sum(injIndex);
nProd = sum(ProdIndex);
W = reservoirP.W;

box = [repmat([2*meter^3/day 500*meter^3/day], nInj, 1);
       repmat([380*barsa 420*barsa], nProd, 1)];
   
lbu = repmat({box(:,1)./uScale},max(reservoirP.simulationControlSteps),1);
ubu = repmat({box(:,2)./uScale},max(reservoirP.simulationControlSteps),1);

nCells = reservoirP.G.cells.num;
lbxS = zeros(nCells,1)./xScale;
ubxS = ones(nCells,1)./xScale;
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);

lbvS = [zeros(numel(vertcat(W(injIndex).cells)),1);
        -inf(numel(vertcat(W(ProdIndex).cells)),1);
        -inf
        ];
ubvS = [inf(numel(vertcat(W(injIndex).cells)),1);
        zeros(numel(vertcat(W(ProdIndex).cells)),1);
        inf
        ];
lbv = repmat({lbvS},totalPredictionSteps,1);
ubv = repmat({ubvS},totalPredictionSteps,1);

u = arrayfun(@(si)schedule2Controls(si,'uScale',uScale),reservoirP.schedule(0~=diff([0;reservoirP.simulationControlSteps]))','UniformOutput',false);


%{
addpath(genpath('../../optimization/testFunctions'));
[~,~,~,simVars,xs,vs,usliced] = simulateSystemSS(u,ss,[]);

[ ei,fi,vi ] = testProfileGradients(xs,u,vs,ss.step,ci,ss.state,'d',1,'pert',1e-7);


%}

plotSol = @(x,u,v,xd) plotSolutionAdjoint(x,u,v,xd,xScale,vScale,uScale,reservoirP);

[u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbxH',lbx,'ubxH',ubx,'lbu',lbu,'ubu',ubu,...
    'tol',1e-6,'lkMax',4,'max_iter',20,'debugLS',false,'saveIt',false,'condense',false,'computeCrossTerm',false,'plotFunc',plotSol,'plot',false);


