function varargout= simulateSystem_R(x,u,v,ss,varargin)

opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'sLeftSeed',[],'guessX',[],'guessV',[],'xRightSeed',[],'uRightSeed',[],'vRightSeed',[],'simVars',[],'eta',0.9);
opt = merge_options(opt, varargin{:});

nR = numel(ss);

gradients = opt.gradients;

xLeftSeed = opt.xLeftSeed;
if isempty(xLeftSeed)
    xLeftSeed = cell(nR,1);
end
vLeftSeed = opt.vLeftSeed;
if isempty(vLeftSeed)
    vLeftSeed = cell(nR,1);
end
sLeftSeed = opt.sLeftSeed;
guessX = opt.guessX;
if isempty(guessX)
    guessX = cell(nR,1);
end
guessV = opt.guessV;
if isempty(guessV)
    guessV = cell(nR,1);
end
xRightSeed = opt.xRightSeed;
if isempty(xRightSeed)
    xRightSeed = cell(nR,1);
end
uRightSeed = opt.uRightSeed;
vRightSeed = opt.vRightSeed;
if isempty(vRightSeed)
    vRightSeed = cell(nR,1);
end
simVars = opt.simVars;
if isempty(simVars)
    simVars = cell(nR,1);
end
withAlgs = true;


xs = cell(nR,1);
vs = cell(nR,1);
J = cell(nR,1);
converged = cell(nR,1);

usliced = cell(nR,1);

for kr = 1:nR
    
    [xs{kr},vs{kr},J{kr},converged{kr},simVars{kr},usliced{kr}] = simulateSystem(x{kr},u,ss{kr},...
        'gradients',gradients,...
        'xLeftSeed',xLeftSeed{kr},...
        'vLeftSeed',vLeftSeed{kr},...
        'guessX',guessX{kr},...
        'guessV',guessV{kr},...
        'xRightSeed',xRightSeed{kr},...
        'uRightSeed',uRightSeed,...
        'simVars',simVars{kr},...
        'withAlgs',withAlgs);
    
end



Jac = [];
if gradients
    if size(xRightSeed{1},1)==0 && size(sLeftSeed,2)==0  % no seeds given
        [s2,JacS] = realization2s(x,u,v,ss,'partials',true,'eta',opt.eta);
        
        Jac.xJx = cellfun(@(Ji)Ji.xJx ,J,'UniformOutput',false);
        Jac.xJu = cellfun(@(Ji)Ji.xJu ,J,'UniformOutput',false);
        Jac.vJx = cellfun(@(Ji)Ji.vJx ,J,'UniformOutput',false);
        Jac.vJu = cellfun(@(Ji)Ji.vJu ,J,'UniformOutput',false);
        Jac.sJv = JacS.Jv;
        Jac.sJx = JacS.Jx;
        Jac.sJu = JacS.Ju;


    
    elseif size(xRightSeed{1},1) ~=0 && size(sLeftSeed,2)==0  % right seeds given

        [s2,JacS] = realization2s(x,u,v,ss,'partials',true,'eta',opt.eta,'vRightSeed',vRightSeed,'xRightSeed',xRightSeed,'uRightSeed',uRightSeed);
        
        Jac.xJ = cellfun(@(Ji)Ji.xJ ,J,'UniformOutput',false);
        Jac.vJ = cellfun(@(Ji)Ji.vJ ,J,'UniformOutput',false);
        Jac.sJ = JacS.J;

    elseif size(xRightSeed{1},1) ==0 && size(sLeftSeed,2)~=0
        
        [s2,JacS] = realization2s(x,u,v,ss,'partials',true,'eta',opt.eta,'leftSeed',sLeftSeed);
        
        Jac.Jx = cellfun(@(Jr,Jvr)cellfun(@plus,Jr.Jx,Jvr,'UniformOutput',false),J',JacS.Jx,'UniformOutput',false);
        Jac.Ju = cellfun(@(Ji)Ji.Ju ,J','UniformOutput',false);
        Jac.Ju = cellfun(@plus,catAndSum(Jac.Ju),JacS.Ju,'UniformOutput',false);
        
        
        
        Jac.Jv = JacS.Jv;
    else
        error('Not allowed to provide rightSeeds and leftSeeds')
    end
else
    [s2] = realization2s(x,u,v,ss,'partials',false,'eta',opt.eta);
end




converged = all(cell2mat(converged));



varargout{1} = xs;
varargout{2} = vs;
varargout{3} = s2;
varargout{4} = Jac;
varargout{5} = converged;
varargout{6} = simVars;
varargout{7} = usliced;





end


function out = catAndSum(M)

dims = cellfun(@(x)size(x,2),M{1});
M = cellfun(@cell2mat,M,'UniformOutput',false);

if issparse(M{1})
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
