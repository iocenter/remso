function [f,Jac] = sepTarget(xs,u,vs,obj,ss,jobSchedule,work2Job,varargin)
% Evaluates a separable target function
%
% SYNOPSIS:
%  [f,Jac] = sepTarget(xs,u,vs,obj,ss,jobSchedule,work2Job)
%  [f,Jac] = sepTarget(xs,u,vs,obj,ss,jobSchedule,work2Job, 'pn', pv, ...)
%
% PARAMETERS:
%
%   xs - Simulated state output.
%
%   vs - Simulated algebraic state output.
%
%   u - control parameters
%
%  obj - separable objective function
%
%  simF - multiple shooting simulator function
%
%  ss - simulator object
%
%  jobSchedule - Object with parallel workers information
%
%  work2Job - Object with parallel workers information
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%  gradients - if true compute the gradient on the line
%
%  leftSeed - vector for vector-Jacobian product.
%
%  (xRightSeeds,uRightSeeds,vRightSeeds) - Vector for the Jacobian-vector
%  product
%
%   usliced - controls for each shooting interval.
%
% RETURNS:   
%
%   f - function value
%
%   Jac - jacobian
%
% SEE ALSO:
%

opt = struct('gradients',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[],'usliced',[]);
opt = merge_options(opt, varargin{:});



totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

nu = numel(u{1});  %% assuming control dimension preserved


xs = distributeVariables( xs,jobSchedule);
if isempty(opt.usliced)
    usliced = cell(totalPredictionSteps,1);
    for k = 1:totalPredictionSteps
        usliced{k} = u{callArroba(ss.ci,{k})};
    end
    usliced = distributeVariables( usliced,jobSchedule);
else
    usliced = distributeVariables( opt.usliced,jobSchedule);
end
if isempty(vs)
    vs = createEmptyCompositeVar(jobSchedule);
else
    vs = distributeVariables( vs,jobSchedule);
end

if ~isempty(opt.xRightSeeds)
    xRightSeed = distributeVariables(opt.xRightSeeds,jobSchedule);
else
    xRightSeed = cell(1,totalPredictionSteps);
end
if isempty(opt.uRightSeeds)
    uRightSeedSliced = createEmptyCompositeVar(jobSchedule);
elseif isa(opt.uRightSeeds,'Composite')
    uRightSeedSliced = opt.uRightSeeds;
else
    uRightSeedSliced = cell(totalPredictionSteps,1);
    for k = 1:totalPredictionSteps
        uRightSeedSliced{k} = opt.uRightSeeds{ss.ci(k)};
    end
    uRightSeedSliced = distributeVariables(uRightSeedSliced,jobSchedule);
end
if isempty(opt.vRightSeeds)
    vRightSeed = createEmptyCompositeVar(jobSchedule);
else
    vRightSeed = distributeVariables(opt.vRightSeeds,jobSchedule);
end

if ~isempty(opt.leftSeed)
    leftSeed = opt.leftSeed;
else
    leftSeed = [];
end



gradientFlag = opt.gradients;


spmd
    nJobsW = numel(work2Job);
    
    fk = cell(nJobsW,1);
    JacStep = cell(nJobsW,1);
    
    for i = 1:nJobsW
        [fk{i},JacStep{i}] = callArroba(obj{i},{xs{i},usliced{i},vs{i}},...
            'partials',gradientFlag,...
            'xRightSeeds',xRightSeed{i},...
            'uRightSeeds',uRightSeedSliced{i},...
            'vRightSeeds',vRightSeed{i},...
            'leftSeed',leftSeed);
    end
    fW = sum(cat(2,fk{:}),2);
    fW = gplus(fW);
end
f =  fW{1};

Jac = [];
if gradientFlag
    if ~isempty(opt.xRightSeeds)
        spmd
            J = zeros(size(fW,1),1);
            for i = 1:nJobsW
                J = J + JacStep{i}.J;
            end
            J = gplus(J);
        end
        Jac.J = J{1};

    else
        
        if ~isempty(opt.leftSeed)
            Jac.Ju = repmat({zeros(size(opt.leftSeed,1),nu)},1,totalControlSteps);
        else
            Jac.Ju = repmat({zeros(size(f,1),nu)},1,totalControlSteps);
        end
        Jac.Jx = cell(1,totalPredictionSteps);
        Jac.Jv = cell(1,totalPredictionSteps);
    
        % TODO: can this be parallelized?
        for w = 1:numel(jobSchedule.work2Job)
            JacStepW = JacStep{w};
            for j = 1:numel(jobSchedule.work2Job{w})
                k = jobSchedule.work2Job{w}(j);
                [i] = ss.ci(k);
          
                Jac.Ju{i} = Jac.Ju{i} + JacStepW{j}.Ju;
                Jac.Jx{k} = JacStepW{j}.Jx;
                Jac.Jv{k} = JacStepW{j}.Jv;
                
            end
        end
    end
end





end

