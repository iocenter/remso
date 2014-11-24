function obj = NPVStepM(wellSol,netSol,schedule,nCells,varargin)
% Compute net present value of a schedule with well solutions
% Inspired on NPVOW
% This function considers only one step, and include the control as variable

opt     = struct('OilPrice',             1.0 , ...
    'WaterProductionCost',  0.1 , ...
    'WaterInjectionCost',   0.1 , ...
    'DiscountFactor',       0.0 , ...
    'ComputePartials',      false,...
    'scale',                1    ,...
    'leftSeed',[],...
    'sign',1);
opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
ri  = opt.WaterInjectionCost  / stb;
d   = opt.DiscountFactor;


% pressure and saturaton vectors just used for place-holding
p  = zeros(nCells, 1);
sW = zeros(nCells, 1);

dt = schedule.step.val;
time = dt + schedule.time;


sol = wellSol;
qWs  = vertcat(sol.qWs);
qOs  = vertcat(sol.qOs);
injInx  = (vertcat(sol.sign) > 0);
status = vertcat(sol.status);

% Remove closed well.
qWs = qWs(status);
qOs = qOs(status);
injInx = injInx(status);
nW  = numel(qWs);
pBHP = zeros(nW, 1); %place-holder
schVal = zeros(nW, 1); %place-holder


if opt.ComputePartials
    [~, ~, qWs, qOs, ~,~,~] = initVariablesADI(p, sW, qWs, qOs, pBHP,netSol,schVal);
end


prodInx = ~injInx;
obj = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
    spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
    +(rw*prodInx - ri*injInx).*qWs );

if ~(size(opt.leftSeed,2)==0)
   obj.jac = cellfun(@(x)opt.leftSeed*x,obj.jac,'UniformOutput',false); 
end
