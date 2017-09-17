
clear

try
    require adjoint mimetic incomp
catch
    mrstModule add adjoint mimetic incomp
end


% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi ad-core ad-props

here = fileparts(mfilename('fullpath'));
if isempty(here)
    here = pwd();
end


% Include REMSO functionalities
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstDerivated')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstDerivated',filesep,'modules,filesep,'adjoint')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'adjoint')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'multipleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remso')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoCrossSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'singleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'utils')));
addpath(genpath(fullfile(here,filesep,'reservoirData')));

[ reservoirP ] = initReservoir( );

gravity off

% whether or not to show output
verbose = false;
verboseLevel = 0;


totalPredictionSteps = numel(reservoirP.schedule);  % MS intervals

% Schedule partition for each control period and for each simulated step
lastControlSteps = findControlFinalSteps( (1:totalPredictionSteps)' );
controlSchedules = reservoirP.schedule;  %
stepSchedules = reservoirP.schedule;
totalControlSteps = totalPredictionSteps;

ci  = arroba(@controlIncidence,2,{(1:totalPredictionSteps)'});  %1:1  1 to 1


qScale = 5*meter^3/day;
pScale = 20*barsa;


xScale = ones(reservoirP.G.cells.num,1);
vScale = qScale*ones(numel(reservoirP.W),1);
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

stepNPV = arroba(@simpleNPVADI,[6,7],{reservoirP.G, reservoirP.S, reservoirP.W, reservoirP.rock, reservoirP.fluid, controls,'scale',1e-5,'sign',-1},true);
vScale = [vScale;1];



step = cell(totalPredictionSteps,1);
for k=1:totalPredictionSteps
    cik = callArroba(ci,{k});
    step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,controls,stepSchedules(k),reservoirP,...
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



nInj = sum(vertcat(reservoirP.W.sign) == 1);
nProd = sum(vertcat(reservoirP.W.sign) == -1);
W = reservoirP.W;

box = [repmat([300*barsa 700*barsa], nInj, 1);
       repmat([100*barsa 200*barsa], nProd, 1)];
   
lbu = repmat({box(:,1)./uScale},totalControlSteps,1);
ubu = repmat({box(:,2)./uScale},totalControlSteps,1);

nCells = reservoirP.G.cells.num;
lbxS = zeros(nCells,1);
ubxS = ones(nCells,1);
lbx = repmat({lbxS},totalPredictionSteps,1);
ubx = repmat({ubxS},totalPredictionSteps,1);


u = arrayfun(@(si)schedule2Controls(si,controls,'uScale',uScale),reservoirP.schedule','UniformOutput',false);

%load optimalu

%{
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'testFunctions')));
[~,~,~,simVars,xs,vs,usliced] = simulateSystemSS(u,ss,[]);

[ ei,fi,vi ] = testProfileGradients(xs,u,vs,ss.step,ci,ss.state);


%}

[u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbxH',lbx,'ubxH',ubx,'lbu',lbu,'ubu',ubu,...
    'tol',1e-6,'lkMax',4,'max_iter',500,'debugLS',false,'saveIt',false);
