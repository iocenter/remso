function [f,Jac] = sepTarget(xs,u,vs,obj,ss,varargin)
% Evaluates a separable target function
%
% SYNOPSIS:
%  [f,Jac] = sepTarget(xs,u,vs,obj,ss)
%  [f,Jac] = sepTarget(xs,u,vs,obj,ss, 'pn', pv, ...)
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

fk = cell(totalPredictionSteps,1);
jacStep = cell(totalPredictionSteps,1);

if isempty(opt.usliced)
    usliced = cell(totalPredictionSteps,1);
    for k = 1:totalPredictionSteps
        usliced{k} = u{callArroba(ss.ci,{k})};
    end
else
    usliced = opt.usliced;
end
if isempty(vs)
    vs = cell(totalPredictionSteps,1);
end

if ~isempty(opt.xRightSeeds)
    xRightSeed = opt.xRightSeeds;
else
    xRightSeed = cell(1,totalPredictionSteps);
end
if isempty(opt.uRightSeeds)
    uRightSeedSliced = cell(totalPredictionSteps,1);
else
    uRightSeedSliced = cell(totalPredictionSteps,1);
    for k = 1:totalPredictionSteps
        uRightSeedSliced{k} = opt.uRightSeeds{callArroba(ss.ci,{k})};
    end
end
if isempty(opt.vRightSeeds)
    vRightSeed = cell(totalPredictionSteps,1);
else
    vRightSeed = opt.vRightSeeds;
end


if ~isempty(opt.leftSeed)
    leftSeed = opt.leftSeed;
else
    leftSeed = [];
end



gradientFlag = opt.gradients;

%parfor k = 1:totalPredictionSteps
for k = 1:totalPredictionSteps
    [fk{k},jacStep{k}]= callArroba(obj{k},{xs{k},usliced{k},vs{k}},...
            'partials',gradientFlag,...
            'xRightSeeds',xRightSeed{k},...
            'uRightSeeds',uRightSeedSliced{k},...
            'vRightSeeds',vRightSeed{k},...
            'leftSeed',leftSeed);
end
f =  sum(cat(2,fk{:}),2);

Jac = [];
if gradientFlag
    if ~isempty(opt.xRightSeeds)
        Jac.J = zeros(size(f,1),1);
        for k = 1:totalPredictionSteps
            Jac.J = Jac.J + jacStep{k}.J;
        end
    else
        
        if ~isempty(opt.leftSeed)
            Jac.Ju = repmat({zeros(size(opt.leftSeed,1),nu)},1,totalControlSteps);
        else
            Jac.Ju = repmat({zeros(size(f,1),nu)},1,totalControlSteps);
        end
        Jac.Jx = cell(1,totalPredictionSteps);
        Jac.Jv = cell(1,totalPredictionSteps);
        for k = 1:totalPredictionSteps
            cik = callArroba(ss.ci,{k});
            Jac.Ju{cik} = Jac.Ju{cik} + jacStep{k}.Ju;
            Jac.Jx{k} = jacStep{k}.Jx;
            Jac.Jv{k} = jacStep{k}.Jv;
        end
    end
    
end





end

