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

opt = struct('gradients',false,'leftSeed',[],'guessV',[],'guessX',[],'simVars',[],'abortNotConvergent',true,'uRightSeeds',[],'printCounter',true);
opt = merge_options(opt, varargin{:});

%% Process inputs & prepare outputs

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

xs = cell(totalPredictionSteps,1);
vs =  cell(totalPredictionSteps,1);
J = cell(totalPredictionSteps,1);
Jo = cell(totalPredictionSteps,1);

if nargin < 3
    target = [];
end
gradientBacward = opt.gradients;
guessV = opt.guessV;
if isempty(guessV)
    guessV = cell(totalPredictionSteps,1);
end
guessX = opt.guessX;
if isempty(guessX)
    guessX = cell(totalPredictionSteps,1);
end

simVars = opt.simVars;
if isempty(simVars);
    simVars = cell(totalPredictionSteps,1);
end

converged = false(totalPredictionSteps,1);

usliced = cell(totalPredictionSteps,1);
cik = arrayfun(@(k)callArroba(ss.ci,{k}),1:totalPredictionSteps);
usliced(1:totalPredictionSteps) = u(cik);

if isfield(ss,'stepClient')
    step = ss.stepClient;
else
    step = ss.step;
end
xStart = ss.state;
if isempty(opt.simVars)
    opt.simVars = cell(totalPredictionSteps,1);
end
fk = cell(totalPredictionSteps,1);
uRightSeeds = opt.uRightSeeds;
if ~isempty(uRightSeeds) && (size(uRightSeeds{1},1)>0)
   gradientForward = opt.gradients;
   gradientBacward = false;
   xSeed = zeros(size(xStart,1),size(uRightSeeds{1},2));
else
   uRightSeeds = cell(totalControlSteps,1);
   gradientForward = false;
   xSeed = [];

end
xsRightSeeds = cell(totalPredictionSteps,1);
vsRightSeeds = cell(totalPredictionSteps,1);



varargout = cell(1,7);


t0 = tic;
k0 = 0;
for k = 1:totalPredictionSteps
    if opt.printCounter
        [t0,k0] = printCounter(1, totalPredictionSteps, k,'Forward Simulation ',t0,k0);
    end
    
    [xs{k},vs{k},J{k},convergence,simVars{k}] = callArroba(step{k},{xStart,usliced{k}},...
        'gradients',gradientForward,...
        'guessX',guessX{k},...
        'guessV',guessV{k},...
        'simVars',simVars{k},...
        'xRightSeeds',xSeed,...
        'uRightSeeds',uRightSeeds{cik(k)});
    
    
    converged(k) = convergence.converged;
    
    if opt.abortNotConvergent && ~convergence.converged
        if opt.printCounter
            [t0,k0] = printCounter(1, totalPredictionSteps, totalPredictionSteps,'Forward Simulation ',t0,k0); % clean counter;
        end
        varargout{3} = converged;
        return
    end
    if gradientForward
        xsRightSeeds{k} = J{k}.xJ;
        vsRightSeeds{k} = J{k}.vJ;
        xSeed = xsRightSeeds{k};
    end
    
    xStart = xs{k};
    
    % take care of run this just once!. If the condition below is true,
    % this will be calculated during the adjoint evaluation
    if ~opt.gradients && ~isempty(target) && iscell(target);
        [fk{k},Jo{k}]= callArroba(target{k},{xs{k},usliced{k},vs{k}},'partials',gradientForward,'xRightSeeds',xsRightSeeds{k},'uRightSeeds',uRightSeeds{cik(k)},'vRightSeeds',vsRightSeeds{k});
    end
    
end

if ~isempty(target) && ~iscell(target); 
	[f,JacObj] = callArroba(target,{xs,u,vs},'gradients',opt.gradients,'leftSeed',opt.leftSeed,'xRightSeeds',xsRightSeeds,'uRightSeeds',uRightSeeds,'vRightSeeds',vsRightSeeds);
    if gradientForward
        gradU = JacObj.J;
    end
elseif ~opt.gradients && ~isempty(target) && iscell(target);
	if gradientForward
        Jo = cellfun(@(J)J.J,Jo,'UniformOutput',false);
        gradU = catAndSum(Jo);
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
if gradientBacward
    if opt.printCounter
        [t0,k0] = printCounter(totalPredictionSteps,1 , totalPredictionSteps,'Backward Simulation',t0,k0);
    end
    lambdaX = cell(1,totalPredictionSteps);
    lambdaV =  cell(1,totalPredictionSteps);
    
    k = totalPredictionSteps;
    
    if iscell(target)
        [fk{k},JacTar]= callArroba(target{k},{xs{k},usliced{k},vs{k}},...
            'partials',opt.gradients,...
            'leftSeed',opt.leftSeed);
        gradU = repmat({zeros(size(JacTar.Ju))},1,totalControlSteps);
   		cik = callArroba(ss.ci,{k});
    	gradU{cik} = JacTar.Ju;
    else
        JacTar.Jv = JacObj.Jv{k};
        JacTar.Jx = JacObj.Jx{k};
		gradU = JacObj.Ju;
    end
    


    lambdaX{k} = -JacTar.Jx;
    lambdaV{k} = -JacTar.Jv;
    
    
    for k = totalPredictionSteps-1:-1:1
        if opt.printCounter
            [t0,k0] = printCounter(totalPredictionSteps,1 , k,'Backward Simulation ',t0,k0);
        end
        
        if iscell(target)
            [fk{k},JacTar]= callArroba(target{k},{xs{k},usliced{k},vs{k}},...
                'partials',opt.gradients,...
                'leftSeed',opt.leftSeed);
			cik = callArroba(ss.ci,{k});
        	gradU{cik} = gradU{cik} + JacTar.Ju;
        else
            JacTar.Jv = JacObj.Jv{k};
            JacTar.Jx = JacObj.Jx{k};
        end
        
        
        [~,~,JacStep,~,simVars{k+1}] = callArroba(step{k+1},{xs{k},usliced{k+1}},...
            'gradients',true,...
            'xLeftSeed',lambdaX{k+1},...
            'vLeftSeed',lambdaV{k+1},...
            'guessX',guessX{k+1},...
            'guessV',guessV{k+1},...
            'simVars',simVars{k+1});

        cikP = callArroba(ss.ci,{k+1});
            
        gradU{cikP} = gradU{cikP} - JacStep.Ju;
        
        
        lambdaX{k} = -JacTar.Jx + JacStep.Jx;
        lambdaV{k} = -JacTar.Jv;
        
        
    end
	%printCounter(1, totalPredictionSteps, totalPredictionSteps, 'BackwardSimSS');

    k = 0;
    [~,~,JacStep,~,simVars{k+1}] = callArroba(step{k+1},{ss.state,usliced{k+1}},...
        'gradients',true,...
        'xLeftSeed',lambdaX{k+1},...
        'vLeftSeed',lambdaV{k+1},...
        'guessX',guessX{k+1},...
        'guessV',guessV{k+1},...
        'simVars',simVars{k+1});

	cikP = callArroba(ss.ci,{k+1});
    
    gradU{cikP} = gradU{cikP} - JacStep.Ju;
    
    
end

if ~isempty(target)
    if iscell(target)
        f =  sum(cat(2,fk{:}),2);
    else
        % f is already defined
    end
else
    f = [];
end
    
varargout{1} = f;

if opt.gradients
    varargout{2} = gradU;
end
varargout{3} = converged;
varargout{4} = simVars;
varargout{5} = xs;
varargout{6} = vs;
varargout{7} = usliced;
varargout{8} = J;







end

function out = catAndSum(M)
if isempty(M)
    out = 0;
elseif any(cellfun(@issparse,M))
    if isrow(M)
        M = M';
    end
    rows= size(M{1},1);
    blocks = numel(M);
    out = sparse( repmat(1:rows,1,blocks),1:rows*blocks,1)*cell2mat(M);
else
    out = sum(cat(3,M{:}),3);    
end

end
