function varargout= simulateSystemSS_R(u,sss,target,varargin)
% Performs a single shooting simulation on each realization


opt = struct('gradients',false,'leftSeed',[],'guessV',[],'guessX',[],'simVars',[],'abortNotConvergent',false);
opt = merge_options(opt, varargin{:});


ss = sss.ss;
jobSchedule = sss.jobSchedule;
fidW = jobSchedule.fidW;
imMaster = jobSchedule.imMaster;

gradients = opt.gradients;
guessV = opt.guessV;
if isempty(guessV)
    guessV = createEmptyDistributedVar(jobSchedule);
end
guessX = opt.guessX;
if isempty(guessX)
    guessX = createEmptyDistributedVar(jobSchedule);
end
simVars = opt.simVars;
if isempty(simVars)
    simVars = createEmptyDistributedVar(jobSchedule);
end
abortNotConvergent = opt.abortNotConvergent;

%spmd
    [o,~,converged,simVars,xs,vs,usliced] = runSS(u,ss,false,[],guessV,guessX,simVars,abortNotConvergent,fidW);
%end

o = bringVariablesMPI(o,jobSchedule);


% TODO: give outputRisk as an input
if imMaster
s2 = outputRisks(o,'eta',sss.eta,'partials',false);
else
s2 = nan;
end

f = [];
fJac = [];
if nargin > 2 && ~isempty(target) && imMaster
    [ f,fJac] = target(s2,u,'gradients',opt.gradients,'leftSeed',opt.leftSeed);
end

g = [];
if gradients
    
    if imMaster
    [s2,JacO] = outputRisks(o,'eta',sss.eta,'partials',true,'leftSeed',fJac.Js);
    else
    s2 = [];JacO.Jo=[];
    end
    
    JacOJo = distributeVariablesMPI(JacO.Jo,jobSchedule);
    
    %spmd
        [~,go,converged,simVars,xs,vs,usliced] = runSS(u,ss,gradients,JacOJo,guessV,guessX,simVars,abortNotConvergent,fidW)  ;
        
        g = catAndSum(go);
        g = gopMPI('+',g,sss.jobSchedule);
    %end
    if imMaster
    uDims = cellfun(@numel,u);
    g = mat2cell(g,size(g,1),uDims);
    
    g = cellfun(@plus,g,fJac.Ju,'UniformOutput',false);
    end
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



function [o,go,converged,simVars,xs,vs,usliced] = runSS(u,ss,gradients,leftSeed,guessV,guessX,simVars,abortNotConvergent,fidW)

nr = numel(ss);
if isempty(fidW)
    printCounter= false;
    printRef = '\b';
    fid = 1;
else
    printCounter= true;
    fid = fidW;
end

if isempty(leftSeed)
    leftSeed = cell(nr,1);
end

o = cell(nr,1);
go = cell(nr,1);
converged=cell(nr,1);
xs = cell(nr,1);
vs = cell(nr,1);
usliced = cell(nr,1);

for r = 1:nr
    if printCounter
        printRef = sprintf('%d/%d',r,nr);
    end
    [o{r},go{r},converged{r},simVars{r},xs{r},vs{r},usliced{r}] = ...
        simulateSystemSS(u,ss{r},ss{r}.outputF,...
        'gradients',gradients,...
        'leftSeed',leftSeed{r},...
        'guessV',guessV{r},...
        'guessX',guessX{r},...
        'simVars',simVars{r},...
        'abortNotConvergent',abortNotConvergent,...
        'printCounter',printCounter,...
        'fid',fid,...
        'printRef',printRef);
end

end

