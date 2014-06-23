function obj = NPVVOStep(wellSols, schedule,nCells, varargin)
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

% pressure and saturaton vectors just used for place-holding
p  = zeros(nCells, 1);
sW = zeros(nCells, 1);
x  = zeros(nCells, 1);


dt = schedule.step.val;
time = dt + schedule.time;

sol = wellSols;
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
schVal = zeros(nW, 1); %place-holder



if opt.ComputePartials
    [~, ~, ~, qWs, qOs, qGs, ~,~] = ...
        initVariablesADI(p, sW, x, qWs, qOs, qGs, pBHP,schVal);
end

prodInx = ~injInx;
obj = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
    spones(ones(1, nW))*( (-ro*prodInx).*qOs +....
    (-rg*prodInx - rig*injInx).*qGs ...
    +(rw*prodInx - riw*injInx).*qWs );

if ~isempty(opt.leftSeed)
   obj.jac = cellfun(@(x)opt.leftSeed*x,obj.jac,'UniformOutput',false); 
end

