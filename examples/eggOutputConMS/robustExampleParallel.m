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
mrstModule add ad-fi ad-core

% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/procedural'));
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
qwScale = 10*meter^3/day;
qoScale = 10*meter^3/day;
objScale = 1/100000;

% to get initial schedule only
[reservoirP,units] = loadEgg('./reservoirData/',0);
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
    
    lbv = cell(numel(work2Job),1);
    ubv = cell(numel(work2Job),1);
    
    nCells = -1;  % in case numel(work2Job) == 0
    xScale = 0;
    nW = 0;
    totalPredictionSteps = 0;
    
    for r=1:numel(work2Job)
        
        % read the realization
        [reservoirP] = loadEgg('./reservoirData/',work2Job(r));
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
       
  
        
        dd = arroba(@drawdown,[1,2],{reservoirP.fluid,'pscale',pScale},true);
        
        nW = numel(W);
        [algFun] = concatenateMrstTargets([dd,stepNPV],false,[nW,1]);
        
        extraAlgScales = ones(nW+1,1); 
                 
        vScale = [wellSolScales;extraAlgScales];  
        
        step = cell(totalPredictionSteps,1);
        
        vDim = numel(vScale);
        lbvr = -inf(vDim*totalPredictionSteps,1);
        ubvr =  inf(vDim*totalPredictionSteps,1);
        lbvr = mat2cell(lbvr,vDim*ones(totalPredictionSteps,1),1);
        ubvr = mat2cell(ubvr,vDim*ones(totalPredictionSteps,1),1);

        simulator = arroba(@mrstSimulationStep,[1,2],{'stop_if_not_converged',false,'force_step',false},true);
        for k = 1:totalPredictionSteps
            cik = callArroba(ci,{k});
            step{k} = arroba(@mrstStep,...
                [1,2],...
                {...
                simulator,...
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
            
            ubvr{k}(end-nW:end-1) = 0;
            
        end
        
        lbv{r} = lbvr;
        ubv{r} = ubvr;
        
        ss{r}.step = step;
        ss{r}.ci = ci;
        ss{r}.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
        ss{r}.outputF = arroba(@lastNV,[1,2,3],{1},true);

    end  % r=1:numel(work2Job)
    
    
	maxState = struct('pressure',inf*barsa,'sW',1);
    minState = struct('pressure',0*barsa,'sW',0);   
    lbxS = setStateValues( minState,'nCells',nCells,'xScale',xScale);
    ubxS = setStateValues( maxState,'nCells',nCells,'xScale',xScale);
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

totalPredictionSteps = numel(schedule.step.val);
selection = true(totalPredictionSteps,1);
obj = @(s,u,varargin)sumSelectionS(s,u,selection,varargin{:});




%%  Bounds for all variables!
maxProdInput = struct('BHP',450*barsa);
minProdInput = struct('BHP',100*barsa);
maxInjInput = struct('RATE',500*meter^3/day);
minInjInput = struct('RATE',0*meter^3/day);


% Control input bounds for all wells!
[ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
    'maxProd',maxProdInput,'minProd',minProdInput,...
    'maxInj',maxInjInput,'minInj',minInjInput,'useScheduleLims',false);
lbu = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales);
ubu = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales);


totalPredictionSteps = numel(schedule.step.val); 
lbs = -inf(totalPredictionSteps,1);
ubs =  inf(totalPredictionSteps,1);


%%  Initialize from previous solution?
u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);

controlWriter = @(u,i) controlWriterMRST(u,i,controlSchedules,cellControlScales,'filename',['./controls/schedule' num2str(i) '.inc'],'units',units);


%Provide the initial simulation as a guess.
[~,~,~,~,xs,vs,~,~] = simulateSystemSS_R(u,sss,[]);


%% call REMSO
[u,x,v,f,xd,M,simVars] = remso(u,sss,obj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,'lbs',lbs,'ubs',ubs,...
    'tol',1e-2,'lkMax',4,'debugLS',true,'max_iter',500,'debugLS',false,'saveIt',true,'controlWriter',controlWriter,'x',xs,'v',vs);

