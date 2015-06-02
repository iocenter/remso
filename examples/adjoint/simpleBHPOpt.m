
clear

try
    require adjoint mimetic incomp
catch
    mrstModule add adjoint mimetic incomp
end


% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi ad-core ad-props

% Include REMSO functionalities
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstDerivated/modules/adjoint'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/adjoint'));
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));

[ reservoirP ] = initReservoir( )

% whether or not to show output
verbose = false;
verboseLevel = 0;



nInj = sum(vertcat(reservoirP.W.sign) == 1);
nProd = sum(vertcat(reservoirP.W.sign) == -1);


qScale = meter^3/day;
pScale = barsa;




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


simulator = @mrstSimulationStep;

wellSol = reservoirP.state.wellSol;

x0 = stateMrst2stateVector(reservoirP.state,'xScale',xScale);
u = schedule2Controls(reservoirP.schedule(1),'uScale',uScale);
schedule = reservoirP.schedule;

[x,v,Jac,convergence,simVars] = mrstStep(x0,u,simulator,wellSol,schedule(1),reservoirP,...
    'gradients',true,...
    'xLeftSeed',[],...
    'vLeftSeed',[],...
    'xRightSeeds',[],...
    'uRightSeeds',[],...
    'guessX',[],...
    'guessV',[],...
    'xScale',xScale,...
    'vScale',vScale,...
    'uScale',uScale,...
    'saveJacobians',true,...
    'simVars',[],...
    'algFun',[]);



f= @(x0,u)mrstStep(x0,u,simulator,wellSol,schedule(1),reservoirP,...
    'gradients',false,...
    'xLeftSeed',[],...
    'vLeftSeed',[],...
    'xRightSeeds',[],...
    'uRightSeeds',[],...
    'guessX',[],...
    'guessV',[],...
    'xScale',xScale,...
    'vScale',vScale,...
    'uScale',uScale,...
    'saveJacobians',true,...
    'simVars',[],...
    'algFun',[]);



%adOpts = admOptions('i', [1,2],'d',[1,2]);
%jac = admDiffFD(f,1,x0,u,adOpts)


xjx = zeros(numel(x),numel(x0));
vjx = zeros(numel(v),numel(x0));
pert = 0.00001;
for k = 1:numel(x0)
    x0P = x0;
    x0N = x0;
    x0P(k) = x0P(k)+pert;
    x0N(k) = x0N(k)-pert;
    
    [xP,vP] = f(x0P,u);
    [xN,vN] = f(x0N,u);

    
    xjx(:,k) = (xP-xN)/(2*pert);
    vjx(:,k) = (vP-vN)/(2*pert);
end
max(max(abs(Jac.xJx - xjx)))
max(max(abs(Jac.vJx - vjx)))




xju = zeros(numel(x),numel(u));
vju = zeros(numel(v),numel(u));
pert = 0.0001;
for k = 1:numel(u)
    uP = u;
    uP(k) = uP(k)+pert;
    
    [xP,vP] = f(x0,uP);
    
    xju(:,k) = (xP-x)/pert;
    vju(:,k) = (vP-v)/pert;
end
max(max(abs(Jac.xJu - xju)))
max(max(abs(Jac.vJu- vju)))



%xjx = zeros(numel(x),numel(u));
%vjx = zeros(numel(v),numel(u));


% 
% [simRes,reports] = runSchedule(reservoirP.state,...
%     reservoirP.G,...
%     reservoirP.S,...
%     reservoirP.W,...
%     reservoirP.rock,...
%     reservoirP.fluid,...
%     reservoirP.schedule,...
%     'Verbose',  verbose , ...
%     'VerboseLevel', verboseLevel);
% 
% 
% 












%
%
% % Run optimization --------------------------------------------------------
% [simRes, schedule, controls, out] = optimizeObjective(G, S, W, rock, ...
%                                         fluid, state, schedule, ...
%                                         controls, objectiveFunction, ...
%                                         'gradTol',       1e-3, ...
%                                         'objChangeTol',  5e-4, ...
%                                         'VerboseLevel', verboseLevel);
%




controls = initControls(reservoirP.schedule);


target = @(simRes,schedule,varargin)  simpleNPVADI(reservoirP.G,...
                                      reservoirP.S,...
                                      reservoirP.W,...
                                      reservoirP.rock,...
                                      reservoirP.fluid,...
                                      simRes,...
                                      schedule,...
                                      controls, varargin{:});






[o,Jaco,convergence,simVars] = targetMrstStep(x0,u,target,simulator,wellSol,schedule(1),reservoirP,...
    'gradients',true,...
    'xRightSeeds',[],...
    'uRightSeeds',[],...
    'guessX',[],...
    'guessV',[],...
    'xScale',xScale,...
    'vScale',vScale,...
    'uScale',uScale,...
    'saveJacobians',true,...
    'simVars',[]);



ft = @(x0,u) targetMrstStep(x0,u,target,simulator,wellSol,schedule(1),reservoirP,...
    'gradients',false,...
    'xRightSeeds',[],...
    'uRightSeeds',[],...
    'guessX',[],...
    'guessV',[],...
    'xScale',xScale,...
    'vScale',vScale,...
    'uScale',uScale,...
    'saveJacobians',true,...
    'simVars',[]);




%adOpts = admOptions('i', [1,2],'d',[1,2]);
%jac = admDiffFD(f,1,x0,u,adOpts)


ojx = zeros(numel(o),numel(x0));
pert = 0.000001;
for k = 1:numel(x0)
    x0P = x0;
    x0P(k) = x0P(k)+pert;
    
    [oP] = ft(x0P,u);
    
    ojx(:,k) = (oP-o)/pert;
end
max(max(abs(Jaco.Jx - ojx)))




oju = zeros(numel(o),numel(u));
pert = 0.0001;
for k = 1:numel(u)
    uP = u;
    uP(k) = uP(k)+pert;
    
    [oP] = ft(x0,uP);
    
    oju(:,k) = (oP-o)/pert;
end
max(max(abs(Jaco.Ju - oju)))




