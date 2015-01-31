function varargout= simulateSystem(x,u,ss,varargin)
% Performs a multiple shooting simulation
%
% SYNOPSIS:
%  [xs,vs,Jac,converged,simVars,usliced] = simulateSystem(x,u,ss)
%  [xs,vs,Jac,converged,simVars,usliced] = simulateSystem(x,u,ss, 'pn', pv, ...)
% PARAMETERS:
%   x - States in the prediction horizon. Used as initial conditions for
%       the shooting intervals
%
%   u - cellarray containing the controls for each control
%       period.
%
%   ss - A simulator structure, containing all the required
%        information on the model.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   gradients - Set to true if gradients must be computed
%
%   xLeftSeed - LHS for vector-Jacobian multiplication related to x
%
%   vLeftSeed - LHS for vector-Jacobian multiplication related to v
%
%   guessX - Simulated xs guess.
%
%   guessV - Simulated vs guess.
%
%   xRightSeed - RHS for Jacobian-vector multiplication related to x
%
%   vRightSeed - RHS for Jacobian-vector multiplication related to v
%
%   uRightSeed - RHS for Jacobian-vector multiplication related to u
%
%   simVars - Provide the simulation result to hot-start the Jacobian
%             calculation
%
%
% RETURNS:
%
%   xs - Simulated state output.
%
%   vs - Simulated algebraic state output.
%
%   Jac - Estimated objective function value.
%
%   converged - True if simulations converged.
%
%   simVars - Simulation variables result.
%
%   usliced - controls for each shooting interval.
%
% SEE ALSO:
%
%
opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'guessX',[],'guessV',[],'xRightSeed',[],'uRightSeed',[],'simVars',[]);
opt = merge_options(opt, varargin{:});

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

nx = numel(ss.state);
nv = ss.nv;
nu = numel(u{1});
withAlgs = nv > 0;


% Depending on what vector are provided for Jacobian-vector and vector
% Jacobian multiplications the gradient calculator will behave differently
Jac = [];
if opt.gradients
    JacStep = cell(1,totalPredictionSteps);
    if isempty(opt.xLeftSeed) && isempty(opt.xRightSeed)
        xJx = cell(totalPredictionSteps,totalPredictionSteps);
        xJu = cell(totalPredictionSteps,totalControlSteps);
        vJx = cell(totalPredictionSteps,totalPredictionSteps);
        vJu = cell(totalPredictionSteps,totalControlSteps);
    elseif ~isempty(opt.xRightSeed) && ~isempty(opt.xLeftSeed)
        error('not implemented')
    elseif ~isempty(opt.xRightSeed)
        xJ = cell(totalPredictionSteps,1);
        vJ = cell(totalPredictionSteps,1);
    elseif ~isempty(opt.xLeftSeed)
        Jx = repmat({zeros(size(opt.xLeftSeed{1},1),nx)},1,totalPredictionSteps);
        Ju = repmat({zeros(size(opt.xLeftSeed{1},1),nu)},1,totalControlSteps);
    end
else
    
end

converged = false(totalPredictionSteps,1);

%%% avoid "not sliced" variables, this additional variables do not represent copy of information!
xs = cell(totalPredictionSteps,1);
vs =  cell(totalPredictionSteps,1);

%% Prepare Inputs to run the simulations
if isempty(opt.guessV)
    guessV = cell(totalPredictionSteps,1);
else
    guessV = opt.guessV;
end
if isempty(opt.guessX)
    guessX = x;
else
    guessX = opt.guessX;
end
if isempty(opt.xLeftSeed)
    xLeftSeed = cell(1,totalPredictionSteps);
else
    xLeftSeed = opt.xLeftSeed;
end
if isempty(opt.vLeftSeed)
    vLeftSeed = cell(1,totalPredictionSteps);
else
    vLeftSeed = opt.vLeftSeed;
end
if isempty(opt.xRightSeed)
    xRightSeed = cell(totalPredictionSteps,1);
    %%% sort right hand sides!  %% Why?? .... because the first initial
    %%% condition is not a variable! the first x affects only on the second
    %%% xs after simulation
else
    xRightSeed = [{zeros(nx,size(opt.xRightSeed{1},2))};opt.xRightSeed(1:totalPredictionSteps-1)];
end
uRightSeedSliced = cell(totalPredictionSteps,1);
if ~isempty(opt.uRightSeed)
    for k = 1:totalPredictionSteps
        uRightSeedSliced{k} = opt.uRightSeed{callArroba(ss.ci,{k})};
    end
end
usliced = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    usliced{k} = u{callArroba(ss.ci,{k})};
end
step = ss.step;
gradientFlag = opt.gradients;
xStart = [ss.state;x(1:totalPredictionSteps-1)];
if isempty(opt.simVars)
    simVars = cell(totalPredictionSteps,1);
else
    simVars = opt.simVars;
end

% Run the simulations in parallel!
%parfor k = 1:totalPredictionSteps
t0 = tic;
k0 = 0;
for k = 1:totalPredictionSteps
    [t0,k0] = printCounter(1, totalPredictionSteps, k,'MS Simulation Sequential',t0,k0);
    
    
    
    [xs{k},vs{k},JacStep{k},convergence,simVars{k}] = step{k}(xStart{k},usliced{k},...
        'gradients',gradientFlag,...
        'xLeftSeed',xLeftSeed{k},...
        'vLeftSeed',vLeftSeed{k},...
        'guessX',guessX{k},...
        'guessV',guessV{k},...
        'xRightSeeds',xRightSeed{k},...
        'uRightSeeds',uRightSeedSliced{k},...
        'simVars',simVars{k});
    
    
    converged(k) = convergence.converged;
    
end


% extract the gradient information!

if opt.gradients
    if isempty(opt.xRightSeed) && isempty(opt.xLeftSeed)
        for k = 1:totalPredictionSteps
            [i] = feval(ss.ci,k);
            xJu{k,i} = JacStep{k}.xJu;
            if withAlgs
                vJu{k,i} = JacStep{k}.vJu;
            end
            %vJx{k,k} = -zeros(nv,nx);   it is assumed zero
            
            if k>1
                xJx{k,k-1} = JacStep{k}.xJx;
                if withAlgs
                    vJx{k,k-1} = JacStep{k}.vJx;
                end
            elseif k == 1
            else
                error('what?')
            end
        end
        Jac.xJu = xJu;
        Jac.xJx = xJx;
        if withAlgs
            Jac.vJu = vJu;
            Jac.vJx = vJx;
        end
        
    elseif ~isempty(opt.xRightSeed)
        for k = 1:totalPredictionSteps
            xJ{k} =  JacStep{k}.xJ;
            if withAlgs
                vJ{k} = JacStep{k}.vJ;
            end
        end
        Jac.xJ = xJ;
        Jac.vJ = vJ;
        
    else %%% if ~isempty(opt.xLeftSeed)
        for k = 1:totalPredictionSteps
            [i] = ss.ci(k);
            Ju{1,i} = Ju{1,i}+JacStep{k}.Ju;
            
            if k>1
                Jx{1,k-1} = JacStep{k}.Jx;
            elseif k == 1
            else
                error('what?')
            end
        end
        Jac.Ju = Ju;
        Jac.Jx = Jx;
        
    end
end

% Check the convergence condition to issue a warning.
if ~all(converged)
    steps = 1:totalPredictionSteps;
    warning(strcat('Simulate System: Steps failed to converge:',num2str(steps(~converged))));
end

varargout{1} = xs;
varargout{2} = vs;

if opt.gradients
    varargout{3} = Jac;
else
    varargout{3} = [];
end
varargout{4} = converged;
varargout{5} = simVars;
varargout{6} = usliced;





end

