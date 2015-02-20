function obj = NPVStepM(forwardStates,schedule,nCells,varargin)
% Compute net present value of a schedule with well solutions
% Inspired in NPVOW
% This function only changes the inputs, and have additional options

opt     = struct('OilPrice',             1.0 , ...
                 'WaterProductionCost',  0.1 , ...
                 'WaterInjectionCost',   0.1 , ...
                 'DiscountFactor',       0.0 , ...
                 'ComputePartials',      false, ...
                 'tStep' ,               [],...
                 'scale',                1    ,...
                 'leftSeed',[],...
                 'sign',1);
opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
ri  = opt.WaterInjectionCost  / stb;
d   = opt.DiscountFactor;


wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);

% pressure and saturaton vectors just used for place-holding
p  = zeros(nCells, 1);
sW = zeros(nCells, 1);

dts   = schedule.step.val;

tSteps = opt.tStep;
if isempty(tSteps) %do all
    time = 0;
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    time = sum(dts(1:(opt.tStep-1)));
    numSteps = 1;
    dts = dts(opt.tStep);
end

obj = cell(1,numSteps);

for step = 1:numSteps
    sol = wellSols{tSteps(step)};
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


    if opt.ComputePartials
        [~, ~, qWs, qOs, ~] = initVariablesADI(p, sW, qWs, qOs, pBHP);
    end

    dt = dts(step);
    time = time + dt;

    prodInx = ~injInx;
    obj{step} = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
        spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
        +(rw*prodInx - ri*injInx).*qWs );
    
    if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
        obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
    end
    
end
