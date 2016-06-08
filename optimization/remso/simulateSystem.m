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
opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'guessX',[],'guessV',[],'xRightSeed',[],'uRightSeed',[],'simVars',[],'withAlgs',false);
opt = merge_options(opt, varargin{:});


withAlgs = opt.withAlgs;

work2Job = ss.work2Job;
jobSchedule = ss.jobSchedule;

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

nx = numel(ss.state);
uDims =  cellfun(@(z)numel(z),u);


% Depending on what vector are provided for Jacobian-vector and vector
% Jacobian multiplications the gradient calculator will behave differently
Jac = [];
if opt.gradients
    if isempty(opt.xLeftSeed) && isempty(opt.xRightSeed)  %% build full Jacobian
        spmd
            nJobsW = numel(work2Job);
            xJx = cell(nJobsW,totalPredictionSteps);
            xJu = cell(nJobsW,totalControlSteps);
            if withAlgs
                vJx = cell(nJobsW,totalPredictionSteps);
                vJu = cell(nJobsW,totalControlSteps);
            end
        end
    elseif ~isempty(opt.xRightSeed) && ~isempty(opt.xLeftSeed)
        error('not implemented')
    elseif ~isempty(opt.xRightSeed)  %% Jacobian-vector product
        spmd
            nJobsW = numel(work2Job);
            xJ = cell(nJobsW,1);
            if withAlgs
                vJ = cell(nJobsW,1);
            end
        end
    elseif ~isempty(opt.xLeftSeed)  %% vector-jacobian product
        Jx = repmat({zeros(size(opt.xLeftSeed{1},1),nx)},1,totalPredictionSteps);
        Ju = arryfun(@(nu)zeros(size(opt.xLeftSeed{1},1),nu),uDims','UniformOutput',false);
    end
else
    
end



%% Prepare Inputs to run the simulations
if isempty(opt.guessV)
    guessV  = createEmptyCompositeVar(jobSchedule);
else
    guessV = distributeVariables(opt.guessV,jobSchedule);
end
if isempty(opt.guessX)
    guessX = distributeVariables(x,jobSchedule);
else
    guessX = distributeVariables(opt.guessX,jobSchedule);
end
if isempty(opt.xLeftSeed)
    xLeftSeed = createEmptyCompositeVar(jobSchedule);
else
    xLeftSeed = distributeVariables(opt.xLeftSeed,jobSchedule);
end
if isempty(opt.vLeftSeed)
    vLeftSeed = createEmptyCompositeVar(jobSchedule);
else
    vLeftSeed = distributeVariables(opt.vLeftSeed,jobSchedule);
end
if isempty(opt.xRightSeed)
    xRightSeed = createEmptyCompositeVar(jobSchedule);
    
    %%% sort right hand sides!  %% Why?? .... because the first initial
    %%% condition is not a variable! the first x affects only on the second
    %%% xs after simulation
elseif isa(opt.xRightSeed,'Composite')
    xRightSeed = opt.xRightSeed;
    spmd
        if labindex<numlabs
            labSend(xRightSeed{end},labindex+1);
        end
        if labindex > 1
            firstCell = labReceive(labindex-1);
        end
        if labindex == 1
            xRightSeed = [{zeros(nx,size(xRightSeed{1},2))};xRightSeed(1:end-1)];
        else
            xRightSeed = [{firstCell};xRightSeed(1:end-1)];
        end
    end
else
    xRightSeed = [{zeros(nx,size(opt.xRightSeed{1},2))};opt.xRightSeed(1:totalPredictionSteps-1)];
    [ xRightSeed ] = distributeVariables( xRightSeed,jobSchedule);
end
if iscell(opt.uRightSeed)
    uRightSeedSliced = cell(totalPredictionSteps,1);
    for k = 1:totalPredictionSteps
        uRightSeedSliced{k} = opt.uRightSeed{callArroba(ss.ci,{k})};
    end
    uRightSeedSliced = distributeVariables( uRightSeedSliced,jobSchedule);
elseif isempty(opt.uRightSeed)
    uRightSeedSliced = createEmptyCompositeVar(jobSchedule);
else
    error('notImplemented')
end
usliced = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
    usliced{k} = u{callArroba(ss.ci,{k})};
end
usliced = distributeVariables( usliced,jobSchedule);
step = ss.step;
gradientFlag = opt.gradients;

% add the dynamical system initial state in the vector of initial
% conditions
if iscell(x)
    xStart = [ss.state;x(1:totalPredictionSteps-1)];
    xStart = distributeVariables( xStart,jobSchedule);
elseif isa(x,'Composite')
    x0 = ss.state;
    spmd
        xStart = x;
        if labindex<numlabs
            labSend(x{end},labindex+1);
        end
        if labindex > 1
            firstCell = labReceive(labindex-1);
        end
        if labindex == 1
            xStart = [{x0};xStart(1:end-1)];
        else
            xStart = [{firstCell};xStart(1:end-1)];
        end
    end
else
    error('notImplemented')
end
if isempty(opt.simVars)
    simVars = createEmptyCompositeVar(jobSchedule);
else
    simVars = distributeVariables(opt.simVars,jobSchedule);
end

% Run the simulations in parallel!
spmd
    %    printCounter(1, totalPredictionSteps, k, 'ForwardSimMS');
    nJobsW = numel(work2Job);
    
    xs = cell(nJobsW,1);
    vs = cell(nJobsW,1);
    JacStep = cell(nJobsW,1);
    convergence = cell(nJobsW,1);
    converged = true;
    
    for i = 1:nJobsW
        
        [xs{i},vs{i},JacStep{i},convergence{i},simVars{i}] = callArroba(step{work2Job(i)},{xStart{i},usliced{i}},...
            'gradients',gradientFlag,...
            'xLeftSeed',xLeftSeed{i},...
            'vLeftSeed',vLeftSeed{i},...
            'guessX',guessX{i},...
            'guessV',guessV{i},...
            'xRightSeeds',xRightSeed{i},...
            'uRightSeeds',uRightSeedSliced{i},...
            'simVars',simVars{i});
        
        
        converged = converged && convergence{i}.converged;
    end
    convergedAll = gop(@and,converged);
end


% extract the gradient information!
ci = ss.ci;
if opt.gradients
    if isempty(opt.xRightSeed) && isempty(opt.xLeftSeed)
        spmd
            for i = 1:nJobsW
                k = work2Job(i);
                cik = callArroba(ci,{k});
                
                xJu{i,cik} = JacStep{i}.xJu;
                if withAlgs
                    vJu{i,cik} = JacStep{i}.vJu;
                end
                %vJx{k,k} = -zeros(nv,nx);   it is assumed zero
                
                if k>1
                    xJx{i,k-1} = JacStep{i}.xJx;
                    if withAlgs
                        vJx{i,k-1} = JacStep{i}.vJx;
                    end
                elseif k == 1
                else
                    error('what?')
                end
                
            end
        end
        Jac.xJu = xJu;
        Jac.xJx = xJx;
        if withAlgs
            Jac.vJu = vJu;
            Jac.vJx = vJx;
        end
        
    elseif ~isempty(opt.xRightSeed)
        spmd
            for i = 1:nJobsW
                xJ{i} = JacStep{i}.xJ;
                if withAlgs
                    vJ{i} = JacStep{i}.vJ;
                end
            end
        end
        Jac.xJ = xJ;
        if withAlgs
            Jac.vJ = vJ;
	    end        
        
    else
        % TODO: see if it is possible to parallelize and if it is worthy
        for w = 1:numel(ss.jobSchedule.work2Job)
            jacStepW = JacStep{w};
            for j = 1:numel(ss.jobSchedule.work2Job{w})
                k = ss.jobSchedule.work2Job{w}(j);
                [i] = callArroba(ss.ci,{k});
                Ju{1,i} = Ju{1,i}+jacStepW{j}.Ju;
                
                if k>1
                    Jx{1,k-1} = jacStepW{j}.Jx;
                elseif k == 1
                else
                    error('what?')
                end
                
            end
        end
        Jac.Ju = Ju;
        Jac.Jx = Jx;
        
    end
end

% Check the convergence condition to issue a warning.
converged = convergedAll{1};
if ~converged
    stringWarn = '';
    for w = numel(ss.jobSchedule.work2Job)
        if ~converged{w}
            convergencei = convergence{w};
            for j = 1:numel(ss.jobSchedule.work2Job{w})
                if ~convergencei{j}.converged
                    stringWarn = [stringWarn,' ',num2str(ss.jobSchedule.work2Job{w}(j))];
                end
            end
        end
    end
    warning(strcat('Simulate System: Steps failed to converge:',stringWarn));
    
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

