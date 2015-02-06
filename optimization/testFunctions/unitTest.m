function [ errorMax ] = unitTest(u1,ss,objStep,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

opt = struct('totalSteps',3,'debug',false,'feasible',false);
opt = merge_options(opt, varargin{:});

initPool();

if isempty(opt.totalSteps)
    opt.totalSteps = getTotalPredictionSteps(ss);
else
    ss.step = ss.step(1:opt.totalSteps);
end
obj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,objStep,ss,varargin{:});

u = repmat({u1},opt.totalSteps,1);
if opt.feasible
    [~,~,~,~,x,v,~] = simulateSystemSS(u,ss);
else
    x = repmat({ss.state},opt.totalSteps,1);
    [~,v] = simulateSystem(x,u,ss);
    v = cellfun(@(z)z + rand(size(z))-0.5,v,'UniformOutput',false);
end
vDims = cellfun(@(z)size(z,1),v);
withAlgs = sum(vDims)>0;


[ errorMax2 ] = testNonlinearGradient(x,u,v,ss,obj,'debug',opt.debug,'withAlgs',withAlgs);

[ errorMax1 ] = testSimStepGradient(ss.state,u1,ss.step{1},'debug',opt.debug);

[ errorMax3 ] = testLiftOptFunc(x,u,ss,@simulateSystem,withAlgs,'testAdjoint',true,'testFwd',true); 

[ errorMax6 ] = testLiftOptReduction(x,u,v,ss,withAlgs);


[ errorMax7 ] = testSingleShooting(u,ss,objStep);







errorMax = max([errorMax1,errorMax2,errorMax3,errorMax6,errorMax7]);

end

