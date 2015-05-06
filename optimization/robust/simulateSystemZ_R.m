function varargout= simulateSystemZ_R(u,x,v,sss,JacTar,simVars)
%
%  Run an adjoint simulation on the Z function in order to compute a
%  gradient with respect to target
%
%


%% Process inputs & prepare outputs

ss = sss.ss;

outputLambda = nargout >1;
lambdaX = [];
lambdaV = [];

if isfield(JacTar,'Js')
    lambdaS = JacTar.Js;
    
    [~,Js] = realization2s(x,u,v,sss,'partials',true,'leftSeed',lambdaS);
    
    if isfield(JacTar,'Jx')
        tJx = JacTar.Jx;
        sJx = Js.Jx;
        
        Jx = sumJacs(tJx,sJx);
        
    else
        Jx = Js.Jx;
    end
    if isfield(JacTar,'Jv')
        tJv = JacTar.Jv;
        sJv = Js.Jv;
        
        Jv = sumJacs(tJv,sJv);
        
    else
        Jv = Js.Jv;
    end
    if isfield(JacTar,'Ju')
        uDims = cellfun(@numel,u);
        JacTar.Ju = sumJacsS(JacTar.Ju,mat2cell(Js.Ju,size(Js.Ju,1),uDims));
    else
        JacTar.Ju = Js.Ju;
    end
else
    lambdaS = [];
end

if outputLambda
    [~,g,~,lambdaX,lambdaV] = runSimulateSystemZ(u,x,v,ss,simVars,Jx,Jv);
else
    [~,g] = runSimulateSystemZ(u,x,v,ss,simVars,Jx,Jv);
end

gradU = catAndSum(g);






if isfield(JacTar,'Ju')
    gradU = gradU + cell2mat(JacTar.Ju);
end
uDims = cellfun(@numel,u);
gradU = mat2cell(gradU,size(gradU,1),uDims);

varargout{1} = gradU;
varargout{2} = lambdaX;
varargout{3} = lambdaV;
varargout{4} = lambdaS;


end



function js = sumJacs(JT,J)

js = cellfun(@sumJacsS,JT,J,'UniformOutput',false);

end
function js = sumJacsS(JT,J)

js = cellfun(@plus,JT,J,'UniformOutput',false);

end

function out = catAndSum(M)
if ~isempty(M)
    if iscell(M{1})
        M = cellfun(@cell2mat,M,'UniformOutput',false);
    end
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

function [f,g,usliced,lambdaX,lambdaV] = runSimulateSystemZ(u,x,v,ss,simVars,Jx,Jv)

nr = numel(ss);

printCounter= true;
fid = 1;
    
JacTarW = cell(nr,1);
f = cell(nr,1);
g = cell(nr,1);
usliced = cell(nr,1);
    
if  nargout >3;
    lambdaX = cell(nr,1);
    lambdaV = cell(nr,1);
    for r = 1:nr
        if printCounter
            printRef = sprintf('%d/%d',r,nr);
        end
        JacTarW{r} = subsasgn(JacTarW{r},struct('type','.','subs','Jx'),Jx{r});
        JacTarW{r} = subsasgn(JacTarW{r},struct('type','.','subs','Jv'),Jv{r});
        [f{r},g{r},usliced{r},lambdaX{r},lambdaV{r}] =...
            simulateSystemZ(u,x{r},v{r},ss{r},[],...
            'simVars',simVars{r},...
            'JacTar',JacTarW{r},...
            'withAlgs',true,...
            'printCounter',printCounter,...
            'fid',fid,...
            'printRef',printRef);
    end
else
    for r = 1:nr
        if printCounter
            printRef = sprintf('%d/%d',r,nr);
        end
        JacTarW{r} = subsasgn(JacTarW{r},struct('type','.','subs','Jx'),Jx{r});
        JacTarW{r} = subsasgn(JacTarW{r},struct('type','.','subs','Jv'),Jv{r});
        [f{r},g{r},usliced{r}] =...
            simulateSystemZ(u,x{r},v{r},ss{r},[],...
            'simVars',simVars{r},...
            'JacTar',JacTarW{r},...
            'withAlgs',true,...
            'printCounter',printCounter,...
            'fid',fid,...
            'printRef',printRef);
    end
    
    lambdaX = [];
    lambdaV = [];
    
end

end
