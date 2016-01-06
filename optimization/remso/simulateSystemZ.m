function varargout= simulateSystemZ(u,x,v,ss,target,varargin)
%
%  Run an adjoint simulation on the Z function in order to compute a
%  gradient with respect to target
%
%
opt = struct('leftSeed',[],'simVars',[],'JacTar',[],'withAlgs',false,'printCounter',true,'fid',1,'printRef','');
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
    [t0,k0] = printCounter(totalPredictionSteps,1 , totalPredictionSteps,['Backward Simulation',opt.printRef,' '],t0,k0,opt.fid);
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
someActive = false;
for k = totalPredictionSteps-1:-1:1
    if opt.printCounter
        [t0,k0] = printCounter(totalPredictionSteps,1 , k,['Backward Simulation',opt.printRef,' '],t0,k0,opt.fid);
    end
    active = any([lambdaX,lambdaV],2);
    if someActive || any(active)
        someActive = true;
        [~,~,JacStep] = callArroba(step{k+1},{x{k},usliced{k+1}},...
            'gradients',true,...
            'xLeftSeed',lambdaX(active,:),...
            'vLeftSeed',lambdaV(active,:),...
            'simVars',simVars{k+1});      
    end
    
    cikP = callArroba(ss.ci,{k+1});
    
    if someActive
        gradU{cikP}(active,:) = gradU{cikP}(active,:) + JacStep.Ju;
    end
    lambdaX = JacTar.Jx{k};
    if someActive
        lambdaX(active,:) = lambdaX(active,:) + JacStep.Jx;
    end
    if isfield(JacTar,'Jv')
        lambdaV = +JacTar.Jv{k};
    end
    if outputLambda
        lambdaXOut{k} = lambdaX;
        lambdaVOut{k} = lambdaV;
    end
    
end

active = any([lambdaX,lambdaV],2);
k = 0;
if someActive || any(active)
    someActive = true;
    [~,~,JacStep] = callArroba(step{k+1},{ss.state,usliced{k+1}},...
        'gradients',true,...
        'xLeftSeed',lambdaX(active,:),...
        'vLeftSeed',lambdaV(active,:),...
        'simVars',simVars{k+1});
end
cikP = callArroba(ss.ci,{k+1});
if someActive
    gradU{cikP}(active,:) = gradU{cikP}(active,:) + JacStep.Ju;
end


varargout = cell(1,6);
varargout{1} = f;
varargout{2} = gradU;
varargout{3} = usliced;
varargout{4} = lambdaXOut;
varargout{5} = lambdaVOut;






end

