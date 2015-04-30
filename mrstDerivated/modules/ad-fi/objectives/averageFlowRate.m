function obj = averageFlowRate(forwardStates,schedule,nCells,varargin)
% Compute net present value of a schedule with well solutions
% Inspired in NPVOW
% This function only changes the inputs, and have additional options

opt     = struct('ComputePartials',false,'leftSeed',[],'scale',1);
opt     = merge_options(opt, varargin{:});


wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);

% pressure and saturaton vectors just used for place-holding
p  = zeros(nCells, 1);
sW = zeros(nCells, 1);

dts   = schedule.step.val;

numSteps = numel(dts);
tSteps = (1:numSteps)';


obj = cell(1,numSteps);

T = sum(dts);
for step = 1:numSteps
    sol = wellSols{tSteps(step)};
    qWs  = vertcat(sol.qWs);
    qOs  = vertcat(sol.qOs);
    injInx  = (vertcat(sol.sign) > 0);
%{
Don't remove closed wells, the gradients size will not be consistent!    
    status = vertcat(sol.status);

    % Remove closed well.
    qWs = qWs(status);
    qOs = qOs(status);
    injInx = injInx(status);
%}
    nW  = numel(qWs);
    pBHP = zeros(nW, 1); %place-holder


    if opt.ComputePartials
        [~, ~, qWs, qOs, ~] = initVariablesADI(p, sW, qWs, qOs, pBHP);
    end

    dt = dts(step);

    prodInx = ~injInx;
    obj{step} = opt.scale*dt/T*[+qOs(injInx)  + qWs(injInx);
                              -qOs(prodInx) - qWs(prodInx)];
    
    if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
        obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
    end
    
end
