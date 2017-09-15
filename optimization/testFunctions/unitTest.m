function [ errorMax,eCross ] = unitTest(u,ss,obj,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

opt = struct('totalSteps',3,'debug',false,'feasible',false,'noise',true);
opt = merge_options(opt, varargin{:});

initPool();

if isempty(opt.totalSteps)
    opt.totalSteps = getTotalPredictionSteps(ss);
else
    ss.step = ss.step(1:opt.totalSteps);
end
if iscell(obj)  %% objective is given as a separable sum
    obj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,varargin{:});
end

% TODO: u may have diferent dimensions, correct!
iMax = callArroba(ss.ci,{opt.totalSteps});
u = u(1:iMax);
if opt.feasible
    [~,~,~,~,x,v,~] = simulateSystemSS(u,ss);

   if opt.noise
       x = cellfun(@(xi)xi+(rand(size(xi))-0.5)*0.01,x,'UniformOutput',false); 
       v = cellfun(@(xi)xi+(rand(size(xi))-0.5)*0.01,v,'UniformOutput',false); 
    end 
else
    x = repmat({ss.state},opt.totalSteps,1);
    [~,v] = simulateSystem(x,u,ss);
    v = cellfun(@(z)z + rand(size(z))-0.5,v,'UniformOutput',false);
    if opt.noise
       x = cellfun(@(xi)xi+(rand(size(xi))-0.5)*0.01,x,'UniformOutput',false); 
       v = cellfun(@(xi)xi+(rand(size(xi))-0.5)*0.01,v,'UniformOutput',false); 
    end 
end
vDims = cellfun(@(z)size(z,1),v);
withAlgs = sum(vDims)>0;


[ errorMax1 ] = testSimStepGradient(ss.state,u{1},ss.step{1},'debug',opt.debug);


[ errorMax2,eCross ] = testNonlinearGradient(x,u,v,ss,obj,'debug',opt.debug,'withAlgs',withAlgs);


[ errorMax3 ] = testLiftOptFunc(x,u,ss,@simulateSystem,withAlgs,'testAdjoint',true,'testFwd',true); 

[ errorMax6 ] = testLiftOptReduction(x,u,v,ss,withAlgs);


[ errorMax7 ] = testSingleShooting(u,ss,obj);







errorMax = max([errorMax1,errorMax2,errorMax3,errorMax6,errorMax7]);

end

