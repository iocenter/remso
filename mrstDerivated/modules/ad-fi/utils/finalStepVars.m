function obj = finalStepVars(step,finalState,wellSol,schedule,finalTime, varargin)
% the final state of the simulation (finalState) should equal stateNext

opt     = struct('ComputePartials',false,'xvScale',[],'xLeftSeed',[],'vLeftSeed',[]);

opt     = merge_options(opt, varargin{:});


p   = finalState.pressure;
sW  = finalState.s(:,1); 
qWs = vertcat(wellSol.qWs);
qOs = vertcat(wellSol.qOs);
pBH = vertcat(wellSol.bhp);


if opt.ComputePartials
        [p, sW, qWs, qOs, pBH] = initVariablesADI(p, sW, qWs, qOs, pBH);
end


dts   = schedule.step.val;
time = sum(dts(1:(step)));
if isfield(schedule,'time')
    time = time + schedule.time;
end

if finalTime ~= time
    p = 0*p;
    sW = 0*sW;
    pBH = 0*pBH;
    qWs = 0*qWs;
    qOs = 0*qOs;
end

obj = [p; sW; qWs; qOs; pBH];

if ~isempty(opt.xvScale)  
   obj = obj./[opt.xvScale];
end

if ~(size(opt.xLeftSeed,2)==0)
   obj.jac = cellfun(@(x)[opt.xLeftSeed,opt.vLeftSeed]*x,obj.jac,'UniformOutput',false); 
end



end




