function [ errorMax ] = unitTest(u1,ss,objStep,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

opt = struct('totalSteps',3,'debug',false);
opt = merge_options(opt, varargin{:});

[status,machineName] = system('hostname');
machineName = machineName(1:end-1);
if (matlabpool('size') == 0)
    if(strcmp(machineName,'tucker') == 1) % checking to see if my pool is already open
        matlabpool open 12
    elseif(strcmp(machineName,'apis') == 1)
%        matlabpool open 8
    elseif(strcmp(machineName,'honoris') == 1)
        matlabpool open 8
    elseif(strcmp(machineName,'beehive') == 1)
        matlabpool open 8
    elseif(strcmp(machineName,'ITK-D1000') == 1)
        matlabpool open 8
	elseif(strcmp(machineName,'dantzig') == 1)
        matlabpool open 12
    else
        machineName = machineName(1:end-1);
        if(strcmp(machineName,'service') == 1)  %% vilje
            matlabpool open 12
        end
    end
end

if isempty(opt.totalSteps)
    opt.totalSteps = getTotalPredictionSteps(ss);
else
    ss.step = ss.step(1:opt.totalSteps);
    obj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,objStep,ss,varargin{:});
end

x = repmat({ss.state},opt.totalSteps,1);
v = repmat({rand(ss.nv,1)},opt.totalSteps,1);
u = repmat({u1},opt.totalSteps,1);

%errorMax2 = 0
[ errorMax2 ] = testNonlinearGradient(x,u,v,ss,obj,'debug',opt.debug);


% deprecated test
%[ errorMax7 ] = testLiftOptFuncAdj(x,u,ss,obj,@simulateSystemResidual);


[ errorMax1 ] = testSimStepGradient(ss.state,u1,ss.step{1},'debug',opt.debug);



[ errorMax3 ] = testLiftOptFunc(x,u,ss,@simulateSystem,'testAdjoint',true,'testFwd',true); 

[ errorMax6 ] = testLiftOptReduction(x,u,v,ss);

[ errorMax4 ] = testLiftOptFunc(x,u,ss,@simulateSystemResidual,'testAdjoint',true,'testFwd',true);


[ errorMax7 ] = testSingleShooting(u,ss,objStep);







errorMax = max([errorMax1,errorMax2,errorMax3,errorMax4,errorMax6,errorMax7]);

end

