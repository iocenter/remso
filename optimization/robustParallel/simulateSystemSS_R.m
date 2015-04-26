function varargout= simulateSystemSS_R(u,sss,target,varargin)
% Performs a single shooting simulation on each realization


opt = struct('gradients',false,'leftSeed',[],'guessV',[],'guessX',[],'simVars',[],'abortNotConvergent',false,'eta',0.9);
opt = merge_options(opt, varargin{:});


ss = sss.ss;
jobSchedule = sss.jobSchedule;

gradients = opt.gradients;
guessV = opt.guessV;
if isempty(guessV)
    guessV = createEmptyCompositeVar(jobSchedule);
end
guessX = opt.guessX;
if isempty(guessX)
    guessX = createEmptyCompositeVar(jobSchedule);
end
simVars = opt.simVars;
if isempty(simVars)
    simVars = createEmptyCompositeVar(jobSchedule);
end
abortNotConvergent = opt.abortNotConvergent;

spmd    
    [o,~,converged,simVars,xs,vs,usliced] = runSS(u,ss,false,[],guessV,guessX,simVars,abortNotConvergent)  ;  
end

o = bringVariables(o,jobSchedule);


% TODO: give outputRisk as an input
s2 = outputRisks(o,'eta',opt.eta,'partials',false);

f = [];
if ~isempty(target)
     [ f,fJac] = target(s2,u,'gradients',opt.gradients,'leftSeed',opt.leftSeed);
end

g = [];
if gradients
    
    [s2,JacO] = outputRisks(o,'eta',opt.eta,'partials',true,'leftSeed',fJac.Js);
    
    JacOJo = distributeVariables(JacO.Jo,jobSchedule);
    
    spmd    
        [~,go,converged,simVars,xs,vs,usliced] = runSS(u,ss,gradients,JacOJo,guessV,guessX,simVars,abortNotConvergent)  ;  
        
        g = catAndSum(go);
        g = gplus(g);
    end
    
    g = g{1};
    
    uDims = cellfun(@numel,u);
    g = mat2cell(g,size(g,1),uDims);
        
	g = cellfun(@plus,g,fJac.Ju,'UniformOutput',false);
    
end


varargout = cell(1,7);
varargout{1} = f;

if opt.gradients
    varargout{2} = g;
end

varargout{3} = converged;
varargout{4} = simVars;
varargout{5} = xs;
varargout{6} = vs;
varargout{7} = s2;
varargout{8} = usliced;



end



function out = catAndSum(M)
if ~isempty(M)
    M = cellfun(@cell2mat,M,'UniformOutput',false);
    
    if any(cellfun(@issparse,M))
        if isrow(M)
            M = M';
        end
        rows= size(M{1},1);
        blocks = numel(M);
        out = sparse( repmat(1:rows,1,blocks),1:rows*blocks,1)*cell2mat(M);
    else
        out = sum(cat(3,M{:}),3);
    end
    
else
    out = 0;
end
end



function [o,go,converged,simVars,xs,vs,usliced] = runSS(u,ss,gradients,leftSeed,guessV,guessX,simVars,abortNotConvergent)

if isempty(leftSeed)
	leftSeed = cell(size(ss));
end

    [o,go,converged,simVars,xs,vs,usliced] = cellfun(...
        @(ssr,leftSeedr,guessVr,guessXr,simVarsr)...
        simulateSystemSS(u,ssr,ssr.outputF,...
        'gradients',gradients,...
        'leftSeed',leftSeedr,...
        'guessV',guessVr,...
        'guessX',guessXr,...
        'simVars',simVarsr,...
        'abortNotConvergent',abortNotConvergent),...
        ss,leftSeed,guessV,guessX,simVars,'UniformOutput',false);
    
end

