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
mrstModule clear
mrstModule add deckformat
mrstModule add mimetic adjoint incomp

% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstDerivated/modules/adjoint'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/adjoint'));
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/parallel'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/remsoCrossSequential'));
addpath(genpath('../../optimization/robust'));
addpath(genpath('../../optimization/robustParallel'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));

mrstVerbose off;

% Open a matlab pool depending on the machine availability
initPool('restart',true);


nR = 100;
%%  Who will do what - Distribute the computational effort!
nWorkers = getNumWorkers();
if nWorkers == 0
    nWorkers = 1;
end
[ jobSchedule ] = divideJobsSequentially(nR ,nWorkers);
jobSchedule.nW = nWorkers;

work2Job = Composite();
for w = 1:nWorkers
    work2Job{w} = jobSchedule.work2Job{w};
end

pScale = 5*barsa;
sScale = 0.01;
qScale = 10*meter^3/day;
objScale = 1/10000000;

[reservoirP,units] = loadIncompressibleEgg('./reservoirData/',0);


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
W = reservoirP.W;
injIndex = vertcat(reservoirP.W.sign) == 1;
ProdIndex = vertcat(reservoirP.W.sign) == -1;

lbvS = [zeros(numel(vertcat(W(injIndex).cells)),1);
    -inf(numel(vertcat(W(ProdIndex).cells)),1);
    -inf
    ];
ubvS = [inf(numel(vertcat(W(injIndex).cells)),1);
    zeros(numel(vertcat(W(ProdIndex).cells)),1);
    inf
    ];

reservoirPBase = reservoirP;

spmd
    
    ss = cell(numel(work2Job),1);
    
    lbv = cell(numel(work2Job),1);
    ubv = cell(numel(work2Job),1);
    
    xScale = 0;
    nW = 0;
    totalPredictionSteps = 0;
    
    for r=1:numel(work2Job)
        
        % read the realization
        [reservoirP] = loadIncompressibleEgg('./reservoirData/',work2Job(r));
        
        
        %% Multiple shooting problem set up
        totalPredictionSteps = numel(reservoirP.schedule);  % MS intervals
        stepSchedules = reservoirP.schedule;
        
        % Piecewise linear control -- mapping the step index to the corresponding
        % control
        ci  = arroba(@controlIncidence,2,{reservoirP.simulationControlSteps});
        
        
        %% Variables Scaling
        xScale = sScale*ones(reservoirP.G.cells.num,1);
        vScale = qScale*ones(numel(vertcat(reservoirP.W.cells)),1);

        
        controls = initControls(reservoirP.schedule);
        
        algFun = arroba(@simpleNPVADI,[6,7],{reservoirP.G, reservoirP.S, reservoirP.W, reservoirP.rock, reservoirP.fluid, controls,'scale',objScale,'sign',-1},true);
        vScale = [vScale;1];
        
        step = cell(totalPredictionSteps,1);
        
        lbvr = repmat({lbvS},totalPredictionSteps,1);
        ubvr = repmat({ubvS},totalPredictionSteps,1);
        
        simulator = @mrstSimulationStep;
        for k = 1:totalPredictionSteps
            cik = callArroba(ci,{k});
            step{k} = arroba(@mrstStep,...
                [1,2],...
                {...
                simulator,...
                reservoirP.state.wellSol,...
                stepSchedules(k),...
                reservoirP,...
                'xScale',...
                xScale,...
                'vScale',...
                vScale,...
                'uScale',...
                uScale,...
                'algFun',algFun,...
                'saveJacobians',false...
                },...
                true);
            
        end
        
        lbv{r} = lbvr;
        ubv{r} = ubvr;
        
        ss{r}.step = step;
        ss{r}.ci = ci;
        ss{r}.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
        ss{r}.outputF = arroba(@lastNV,[1,2,3],{1},true);

    end  % r=1:numel(work2Job)
    
    
    nCells = reservoirP.G.cells.num;
    lbxS = zeros(nCells,1)./xScale;
    ubxS = ones(nCells,1)./xScale;
    lbx = repmat({lbxS},totalPredictionSteps,1);
    ubx = repmat({ubxS},totalPredictionSteps,1);
    
    
    outputName = sprintf('w%d.log', labindex);
    fidW = fopen(outputName,'w');
end
jobSchedule.fidW = fidW;
sss.ss = ss;
sss.nR = nR;
sss.jobSchedule = jobSchedule;
sss.eta = 0;


reservoirP= reservoirPBase ;

% to get initial schedule only
schedule = reservoirP.schedule;

injIndex = vertcat(reservoirP.W.sign) == 1;
ProdIndex = vertcat(reservoirP.W.sign) == -1;

totalPredictionSteps = numel(schedule);
selection = true(totalPredictionSteps,1);
obj = @(s,u,varargin)sumSelectionS(s,u,selection,varargin{:});


nInj = sum(injIndex);
nProd = sum(ProdIndex);
W = reservoirP.W;

box = [repmat([0*meter^3/day 500*meter^3/day], nInj, 1);
       repmat([100*barsa 450*barsa], nProd, 1)];
   
lbu = repmat({box(:,1)./uScale},max(reservoirP.simulationControlSteps),1);
ubu = repmat({box(:,2)./uScale},max(reservoirP.simulationControlSteps),1);


lbs = -inf(totalPredictionSteps,1);
ubs =  inf(totalPredictionSteps,1);


%%  Initialize from previous solution?
u = arrayfun(@(si)schedule2Controls(si,'uScale',uScale),reservoirP.schedule(0~=diff([0;reservoirP.simulationControlSteps]))','UniformOutput',false);


loadPrevious = exist('./iterates/itVars_r1.mat','file') ~= 0;
work2Job = jobSchedule.work2Job;
if loadPrevious
	spmd
    nRw = numel(work2Job{labindex});
    
    x = cell(nRw,1);
    xs = cell(nRw,1);
    v  = cell(nRw,1);
    vs = cell(nRw,1);
    simVars = cell(nRw,1);
    for r = 1:numel(work2Job{labindex})
        [u,x{r},xs{r},v{r},vs{r},simVars{r}] = loadItVars('dir','./iterates/','it',0,'r',work2Job{labindex}(r));
    end
    end
    u = u{1};
else
    %Provide the initial simulation as a guess.
	[~,~,~,simVars,xs,vs,~,~] = simulateSystemSS_R(u,sss,[]);
end





if ~loadPrevious
spmd
for r = 1:numel(work2Job{labindex})
    saveItVars(u,xs{r},xs{r},vs{r},vs{r},simVars{r},...
        'dir','./iterates/',...
        'it',0,...
        'r',work2Job{labindex}(r),...
        'keepPreviousIt',true);
end
end
end



%% call REMSO
[u,x,v,f,xd,M,simVars] = remso(u,sss,obj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,'lbs',lbs,'ubs',ubs,...
    'tol',1e-2,'lkMax',4,'debugLS',true,'max_iter',500,'debugLS',false,'saveIt',true,'x',xs,'v',vs,'computeCrossTerm',false);

