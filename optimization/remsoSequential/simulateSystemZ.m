function varargout= simulateSystemZ(u,xd,vd,ss,target,varargin)
%
%  Run an adjoint simulation on the Z function in order to compute a
%  gradient with respect to target
%
%
opt = struct('gradients',false,'leftSeed',[],'guessV',[],'guessX',[],'simVars',[],'abortNotConvergent',true,'JacTar',[]);
opt = merge_options(opt, varargin{:});

%% Process inputs & prepare outputs

totalPredictionSteps = getTotalPredictionSteps(ss);

xs = cell(totalPredictionSteps,1);
vs = cell(totalPredictionSteps,1);
zxs = cell(totalPredictionSteps,1);
zvs =  cell(totalPredictionSteps,1);

withAlgs = ss.nv > 0;
if isempty(opt.guessV)
    opt.guessV = cell(totalPredictionSteps,1);
end
if isempty(opt.guessX)
    opt.guessX = cell(totalPredictionSteps,1);
end

simVarsOut = cell(totalPredictionSteps,1);

converged = false(totalPredictionSteps,1);

usliced = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    usliced{k} = u{callArroba(ss.ci,{k})};
end
step = ss.step;
xStart = ss.state;
if isempty(opt.simVars)
    opt.simVars = cell(totalPredictionSteps,1);
end

varargout = cell(1,7);


t0 = tic;
k0 = 0;
for k = 1:totalPredictionSteps
    [t0,k0] = printCounter(1, totalPredictionSteps, k,'Forward Simulation ',t0,k0);
    
    [xs{k},vs{k},~,convergence,simVarsOut{k}] = step{k}(xStart,usliced{k},...
        'gradients',false,...
        'guessX',opt.guessX{k},...
        'guessV',opt.guessV{k},...
        'simVars',opt.simVars{k});
    
    
    converged(k) = convergence.converged;
    
    if opt.abortNotConvergent && ~convergence.converged
        varargout{3} = false;
        return
    end
    
    zxs{k} = xs{k} - xd{k};
    if withAlgs
        zvs{k} = vs{k} - vd{k};
    end
    xStart = zxs{k};
    
    
end

% take care of run this just once!. If the condition below is true,
% this will be calculated during the adjoint evaluation
if ~isempty(target);
    [f,JacTar] = callArroba(target,{zxs,u,zvs},'gradients',opt.gradients,'usliced',usliced,'leftSeed',opt.leftSeed);
end
if ~isempty(opt.JacTar);
    f = [];
    JacTar = opt.JacTar;
end


% check convergence
if ~all(converged)
    steps = 1:totalPredictionSteps;
    warning(strcat('Simulate System: Steps failed to converge:',num2str(steps(~converged))));
end

% Run the adjoint simulation to get the gradients of the target function!
t0 = tic;
k0 = totalPredictionSteps+1;
if opt.gradients
    [t0,k0] = printCounter(totalPredictionSteps,1 , totalPredictionSteps,'Backward Simulation',t0,k0);
    
    lambdaX = cell(1,totalPredictionSteps);
    lambdaV =  cell(1,totalPredictionSteps);
    
    k = totalPredictionSteps;
    
%     
%     [fk{k},JacTar]= callArroba(target{k},{zxs{k},usliced{k},zvs{k}},...
%         'partials',opt.gradients,...
%         'leftSeed',opt.leftSeed);
%     
    gradU = cellfun(@(xx)zeros(size(xx)),JacTar.Ju,'uniformOutput',false);
    
    lambdaX{k} = +JacTar.Jx{k};
    if isfield(JacTar,'Jv')
        lambdaV{k} = +JacTar.Jv{k};
    end    
    for k = totalPredictionSteps-1:-1:1
        [t0,k0] = printCounter(totalPredictionSteps,1 , k,'Backward Simulation ',t0,k0);
        
%         [fk{k},JacTar]= callArroba(target{k},{zxs{k},usliced{k},zvs{k}},...
%             'partials',opt.gradients,...
%             'leftSeed',opt.leftSeed);
%         
        
        [~,~,JacStep,~,simVarsOut{k+1}] = step{k+1}(zxs{k},usliced{k+1},...
            'gradients',true,...
            'xLeftSeed',lambdaX{k+1},...
            'vLeftSeed',lambdaV{k+1},...
            'guessX',opt.guessX{k+1},...
            'guessV',opt.guessV{k+1},...
            'simVars',simVarsOut{k+1});
        
        cik = callArroba(ss.ci,{k});
        cikP = callArroba(ss.ci,{k+1});
        
        gradU{cikP} = gradU{cikP} + JacStep.Ju;
                
        lambdaX{k} = +JacTar.Jx{k} + JacStep.Jx;
        if isfield(JacTar,'Jv')
            lambdaV{k} = +JacTar.Jv{k};
        end
        
        
    end
    %printCounter(1, totalPredictionSteps, totalPredictionSteps, 'BackwardSimSS');
    
    k = 0;
    [~,~,JacStep,~,simVarsOut{k+1}] = step{k+1}(ss.state,usliced{k+1},...
        'gradients',true,...
        'xLeftSeed',lambdaX{k+1},...
        'vLeftSeed',lambdaV{k+1},...
        'guessX',opt.guessX{k+1},...
        'guessV',opt.guessV{k+1},...
        'simVars',simVarsOut{k+1});
    
    cikP = callArroba(ss.ci,{k+1});
    
    gradU{cikP} = gradU{cikP} + JacStep.Ju;
    
    
end

gradU = cellfun(@plus,gradU,JacTar.Ju,'UniformOutput',false);


varargout{1} = f;

if opt.gradients
    varargout{2} = gradU;
end
varargout{3} = converged;
varargout{4} = simVarsOut;
varargout{5} = xs;
varargout{6} = vs;
varargout{7} = usliced;
varargout{8} = lambdaX;
varargout{9} = lambdaV;






end

