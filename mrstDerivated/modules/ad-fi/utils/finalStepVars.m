function obj = finalStepVars(step,finalState,wellSol,schedule,finalTime, varargin)
% the final state of the simulation (finalState) should equal stateNext

opt     = struct('ComputePartials',false,'xvScale',[],'xLeftSeed',[],'vLeftSeed',[]);

opt     = merge_options(opt, varargin{:});


p   = finalState.pressure;
sW  = finalState.s(:,1); 
qWs = vertcat(wellSol.qWs);
qOs = vertcat(wellSol.qOs);
if isfield(wellSol,'pressure')
	pBH = vertcat(wellSol.pressure);    
else
	pBH = vertcat(wellSol.bhp);
end

if opt.ComputePartials
        [p, sW, qWs, qOs, pBH] = initVariablesADI(p, sW, qWs, qOs, pBH);
end


dts   = schedule.step.val;
time = sum(dts(1:(step)));
if isfield(schedule,'time')
    time = time + schedule.time;
end

if finalTime ~= time
    p = zeros(numel(finalState.pressure))*p;
    sW = zeros(numel(finalState.s(:,1)))*sW;
    pBH = zeros(numel(double(pBH)))*pBH;
    qWs = zeros(numel(double(qWs)))*qWs;
    qOs = zeros(numel(double(qOs)))*qOs;
end

obj = [p; sW; qWs; qOs; pBH];

if ~isempty(opt.xvScale)  
   obj = obj./[opt.xvScale];
end

if ~isempty(opt.xLeftSeed)
   obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false); 
end



end




