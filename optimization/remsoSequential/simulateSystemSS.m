function varargout= simulateSystemSS(u,ss,target,varargin)
% Performs a single shooting simulation
%
% SYNOPSIS:
%  [f,gradU,converged,simVarsOut,xs,vs,usliced] = simulateSystemSS(u,ss,target)
%  [f,gradU,converged,simVarsOut,xs,vs,usliced] = simulateSystemSS(u,ss,target, 'pn', pv, ...)
% PARAMETERS:
%   u - cellarray containing the controls for each control
%       period.
%
%   ss - A simulator structure, containing all the required
%        information on the model.
%
%   target - Separable function of the state, algebraic states and
%            controls. This function is evaluated at each point in the
%            prediction horizon, and the sum of all evaluations is
%            returned.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%           
%   leftSeed - LHS for vector-Jacobian multiplication related to the target
%              function
%
%   guessX - Simulated xs guess.
%
%   guessV - Simulated vs guess.
%
%
%   simVars - Provides the simulation internal data for hotstarting future
%             calls
%
%
% RETURNS:
%
%   f - Sum of the target function evaluations on each predicted step
%   
%   gradU - gradient of the target function w.r.t. u
%
%   converged - true if the simulation converged for all steps
%   
%   simVarsOut - Simulator variables at each simulated step.
%
%   xs - predicted states
%
%   vs - predicted algebraic states
%
%   usliced - controls for each shooting interval.
%
% SEE ALSO:
%
%

opt = struct('gradients',false,'leftSeed',[],'guessV',[],'guessX',[],'simVars',[],'abortNotConvergent',true);
opt = merge_options(opt, varargin{:});

%% Process inputs & prepare outputs

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

xs = cell(totalPredictionSteps,1);
vs =  cell(totalPredictionSteps,1);

if nargin < 3
    target = [];
end

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
fk = cell(totalPredictionSteps,1);

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
    
    
    xStart = xs{k};
    
    % take care of run this just once!. If the condition below is true,
    % this will be calculated during the adjoint evaluation
    if ~opt.gradients && ~isempty(target);
        [fk{k}]= callArroba(target{k},{xs{k},usliced{k},vs{k}},'partials',false);
    end
    
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
    
    
    [fk{k},JacTar]= callArroba(target{k},{xs{k},usliced{k},vs{k}},...
        'partials',opt.gradients,...
        'leftSeed',opt.leftSeed);
    
	gradU = repmat({zeros(size(JacTar.Ju))},1,totalControlSteps);

    lambdaX{k} = -JacTar.Jx;
    lambdaV{k} = -JacTar.Jv;
    
    cik = callArroba(ss.ci,{k});
    gradU{cik} = JacTar.Ju;
    
    for k = totalPredictionSteps-1:-1:1
        [t0,k0] = printCounter(totalPredictionSteps,1 , k,'Backward Simulation ',t0,k0);

        [fk{k},JacTar]= callArroba(target{k},{xs{k},usliced{k},vs{k}},...
            'partials',opt.gradients,...
            'leftSeed',opt.leftSeed);
        
        
        [~,~,JacStep,~,simVarsOut{k+1}] = step{k+1}(xs{k},usliced{k+1},...
            'gradients',true,...
            'xLeftSeed',lambdaX{k+1},...
            'vLeftSeed',lambdaV{k+1},...
            'guessX',opt.guessX{k+1},...
            'guessV',opt.guessV{k+1},...
            'simVars',simVarsOut{k+1});

        cik = callArroba(ss.ci,{k});
        cikP = callArroba(ss.ci,{k+1});
            
        gradU{cikP} = gradU{cikP} - JacStep.Ju;
        gradU{cik} = gradU{cik} + JacTar.Ju;
        
        
        lambdaX{k} = -JacTar.Jx + JacStep.Jx;
        lambdaV{k} = -JacTar.Jv;
        
        
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
    
    gradU{cikP} = gradU{cikP} - JacStep.Ju;
    
    
end

f =  sum(cat(2,fk{:}),2);

varargout{1} = f;

if opt.gradients
    varargout{2} = gradU;
end
varargout{3} = converged;
varargout{4} = simVarsOut;
varargout{5} = xs;
varargout{6} = vs;
varargout{7} = usliced;







end

