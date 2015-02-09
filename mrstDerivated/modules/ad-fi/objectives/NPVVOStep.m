function obj = NPVVOStep(forwardStates,schedule,nCells,varargin)
% Compute net present value of a schedule with well solutions
% Inspired on NPVVO
% This function considers only one step, and include the control as variable


opt     = struct('OilPrice',             1.0 , ...
    'GasPrice',             0.1 , ...
    'GasInjectionCost',     0.1 , ...
    'WaterProductionCost',  0.1 , ...
    'WaterInjectionCost',   0.1 , ...
    'DiscountFactor',       0.0 , ...
    'ComputePartials',      false, ...
    'scale',                1    ,...
    'leftSeed',[],...
    'sign',1);
opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
riw  = opt.WaterInjectionCost / stb;
rg  = opt.GasPrice   / stb;
rig  = opt.GasInjectionCost   / stb;

d   = opt.DiscountFactor;

wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);

% pressure and saturaton vectors just used for place-holding
p  = zeros(nCells, 1);
sW = zeros(nCells, 1);
x  = zeros(nCells, 1);

dts = schedule.step.val;
time = 0;
numSteps = numel(dts);
tSteps = (1:numSteps)';


obj = cell(1,numSteps);

for step = 1:numSteps
    sol = wellSols{tSteps(step)};
    qWs  = vertcat(sol.qWs);
    qOs  = vertcat(sol.qOs);
    qGs  = vertcat(sol.qGs);
    injInx  = (vertcat(sol.sign) > 0);
    status = vertcat(sol.status);
    
    % Remove closed well.
    qWs = qWs(status);
    qOs = qOs(status);
    qGs = qGs(status);
    injInx = injInx(status);
    nW  = numel(qWs);
    pBHP = zeros(nW, 1); %place-holder
    
    
    
    if opt.ComputePartials
        [~, ~, ~, qWs, qOs, qGs, ~] = ...
            initVariablesADI(p, sW, x, qWs, qOs, qGs, pBHP);
    end
    
    dt = dts(step);
    time = time + dt;
    prodInx = ~injInx;
    obj{step} = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
        spones(ones(1, nW))*( (-ro*prodInx).*qOs +....
        (-rg*prodInx - rig*injInx).*qGs ...
        +(rw*prodInx - riw*injInx).*qWs );
    
    if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
        obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
    end
    
    
end