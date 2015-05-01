function varargout= simulateSystemZ(u,x,v,ss,target,varargin)
%
%  Run an adjoint simulation on the Z function in order to compute a
%  gradient with respect to target
%
%
opt = struct('leftSeed',[],'simVars',[],'JacTar',[],'withAlgs',false,'printCounter',true);
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

outputLambda = nargout > 3;
lambdaXOut = cell(1,totalPredictionSteps);
lambdaVOut =  cell(1,totalPredictionSteps);



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


if isfield(JacTar,'Ju') && ~isempty(JacTar.Ju) 
    gradU = JacTar.Ju;
else
    % the target is independent of u
    uDims = cellfun(@numel,u);
    gradU = mat2cell(zeros(size(JacTar.Jx{k},1),sum(uDims)),size(JacTar.Jx{k},1),uDims);    
end

if size(JacTar.Jx{k},1) ==  0   %% deal with this special case that might be often in robust opt
	varargout = cell(1,6);
	varargout{1} = f;
	varargout{2} = gradU;
    varargout{3} = usliced;
    varargout{4} = JacTar.Jx;  %% these are empty and with the right size!
    varargout{5} = JacTar.Jv;  %% these are empty and with the right size!
    return
end


% Run the adjoint simulation to get the gradients of the target function!
if opt.printCounter
    t0 = tic;
    k0 = totalPredictionSteps+1;
    [t0,k0] = printCounter(totalPredictionSteps,1 , totalPredictionSteps,'Backward Simulation',t0,k0);
end


k = totalPredictionSteps;

lambdaX = +JacTar.Jx{k};
if isfield(JacTar,'Jv')
    lambdaV = +JacTar.Jv{k};
end
if outputLambda
	lambdaXOut{k} = lambdaX;
    lambdaVOut{k} = lambdaV;
end

for k = totalPredictionSteps-1:-1:1
    if opt.printCounter
        [t0,k0] = printCounter(totalPredictionSteps,1 , k,'Backward Simulation ',t0,k0);
    end
    
    [~,~,JacStep] = callArroba(step{k+1},{x{k},usliced{k+1}},...
        'gradients',true,...
        'xLeftSeed',lambdaX,...
        'vLeftSeed',lambdaV,...
        'simVars',simVars{k+1});
    
    cikP = callArroba(ss.ci,{k+1});
    
    gradU{cikP} = gradU{cikP} + JacStep.Ju;
    
    lambdaX = +JacTar.Jx{k} + JacStep.Jx;
    if isfield(JacTar,'Jv')
        lambdaV = +JacTar.Jv{k};
    end
    if outputLambda
        lambdaXOut{k} = lambdaX;
        lambdaVOut{k} = lambdaV;
    end    
    
end

k = 0;
[~,~,JacStep] = callArroba(step{k+1},{ss.state,usliced{k+1}},...
    'gradients',true,...
    'xLeftSeed',lambdaX,...
    'vLeftSeed',lambdaV,...
    'simVars',simVars{k+1});

cikP = callArroba(ss.ci,{k+1});
gradU{cikP} = gradU{cikP} + JacStep.Ju;



varargout = cell(1,6);
varargout{1} = f;
varargout{2} = gradU;
varargout{3} = usliced;
varargout{4} = lambdaXOut;
varargout{5} = lambdaVOut;






end

