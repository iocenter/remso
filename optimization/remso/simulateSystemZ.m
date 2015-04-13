function varargout= simulateSystemZ(u,x,v,ss,target,varargin)
%
%  Run an adjoint simulation on the Z function in order to compute a
%  gradient with respect to target
%
%
opt = struct('leftSeed',[],'simVars',[],'JacTar',[],'withAlgs',false);
opt = merge_options(opt, varargin{:});

%% Process inputs & prepare outputs

totalPredictionSteps = getTotalPredictionSteps(ss);

withAlgs = opt.withAlgs;

usliced = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    usliced{k} = u{callArroba(ss.ci,{k})};
end
if isfield(ss,'stepClient')
    step = ss.stepClient;
else
    step = ss.step;
end
if isempty(opt.simVars)
    simVars = cell(totalPredictionSteps,1);
elseif iscell(opt.simVars)
    simVars = opt.simVars;
else
    simVars = bringVariables(opt.simVars,ss.jobSchedule);
end




% take care of run this just once!. If the condition below is true,
% this will be calculated during the adjoint evaluation
if ~isempty(opt.JacTar);
    f = [];
    JacTar = opt.JacTar;
else 
    if withAlgs
        [f,JacTar] = callArroba(target,{x,u,v},'gradients',true,'usliced',usliced,'leftSeed',opt.leftSeed);
    else
        [f,JacTar] = callArroba(target,{x,u},'gradients',true,'usliced',usliced,'leftSeed',opt.leftSeed);
    end
end



gradU = JacTar.Ju;


% Run the adjoint simulation to get the gradients of the target function!

t0 = tic;
k0 = totalPredictionSteps+1;
[t0,k0] = printCounter(totalPredictionSteps,1 , totalPredictionSteps,'Backward Simulation',t0,k0);

lambdaX = cell(1,totalPredictionSteps);
lambdaV =  cell(1,totalPredictionSteps);

k = totalPredictionSteps;

lambdaX{k} = +JacTar.Jx{k};
if isfield(JacTar,'Jv')
    lambdaV{k} = +JacTar.Jv{k};
end


for k = totalPredictionSteps-1:-1:1
    [t0,k0] = printCounter(totalPredictionSteps,1 , k,'Backward Simulation ',t0,k0);
    
    [~,~,JacStep,~,simVars{k+1}] = step{k+1}(x{k},usliced{k+1},...
        'gradients',true,...
        'xLeftSeed',lambdaX{k+1},...
        'vLeftSeed',lambdaV{k+1},...
        'simVars',simVars{k+1});
    
    cikP = callArroba(ss.ci,{k+1});
    
    gradU{cikP} = gradU{cikP} + JacStep.Ju;
    
    lambdaX{k} = +JacTar.Jx{k} + JacStep.Jx;
    if isfield(JacTar,'Jv')
        lambdaV{k} = +JacTar.Jv{k};
    end
    
    
end

k = 0;
[~,~,JacStep,~,simVars{k+1}] = step{k+1}(ss.state,usliced{k+1},...
    'gradients',true,...
    'xLeftSeed',lambdaX{k+1},...
    'vLeftSeed',lambdaV{k+1},...
    'simVars',simVars{k+1});

cikP = callArroba(ss.ci,{k+1});
gradU{cikP} = gradU{cikP} + JacStep.Ju;



varargout = cell(1,6);
varargout{1} = f;
varargout{2} = gradU;
varargout{3} = simVars;
varargout{4} = usliced;
varargout{5} = lambdaX;
varargout{6} = lambdaV;






end

