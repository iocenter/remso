function varargout= simulateSystemSS_R(u,ss,target,varargin)
% Performs a single shooting simulation on each realization


opt = struct('gradients',false,'leftSeed',[],'guessV',[],'guessX',[],'simVars',[],'abortNotConvergent',false,'eta',0.9);
opt = merge_options(opt, varargin{:});


nR = numel(ss);


gradients = opt.gradients;
leftSeed = opt.leftSeed;
if isempty(leftSeed)
    leftSeed = cell(nR,1);
end
guessV = opt.guessV;
if isempty(guessV)
    guessV = cell(nR,1);
end
guessX = opt.guessX;
if isempty(guessX)
    guessX = cell(nR,1);
end
simVars = opt.simVars;
if isempty(simVars)
    simVars = cell(nR,1);
end
abortNotConvergent = opt.abortNotConvergent;


o = cell(nR,1);
go = cell(nR,1);
f = [];
g = [];
converged = cell(nR,1);
xs = cell(nR,1);
vs =  cell(nR,1);
usliced = cell(nR,1);

for kr = 1:nR
    
    
    [o{kr},~,converged{kr},simVars{kr},xs{kr},vs{kr},usliced{kr}] = simulateSystemSS(u,ss{kr},ss{kr}.outputF,...
        'gradients',false,...
        'leftSeed',leftSeed{kr},...
        'guessV',guessV{kr},...
        'guessX',guessX{kr},...
        'simVars',simVars{kr},...
        'abortNotConvergent',abortNotConvergent);
    
end

% TODO: give outputRisk as an input
s2 = outputRisks(o,'eta',opt.eta,'partials',false);

if ~isempty(target)
     [ f,fJac] = target(s2,u,'gradients',opt.gradients,'leftSeed',opt.leftSeed);
end

if opt.gradients
    
    [s2,JacO] = outputRisks(o,'eta',opt.eta,'partials',true,'leftSeed',fJac.Js);
    
     for kr = 1:nR
        [~,go{kr},converged{kr},simVars{kr},xs{kr},vs{kr},usliced{kr}] = simulateSystemSS(u,ss{kr},ss{kr}.outputF,...
            'gradients',true,...
            'leftSeed',JacO.Jo{kr},...
            'guessV',guessV{kr},...
            'guessX',guessX{kr},...
            'simVars',simVars{kr},...
            'abortNotConvergent',abortNotConvergent);  
     end
    
        g = catAndSum(go);
        
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

dims = cellfun(@(x)size(x,2),M{1});
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

out = mat2cell(out,size(out,1),dims);

end
