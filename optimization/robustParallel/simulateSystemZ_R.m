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
        spmd
        Jx = sumJacs(tJx,sJx);
        end
    else
        Jx = Js.Jx;
    end
    if isfield(JacTar,'Jv')
        tJv = JacTar.Jv;
        sJv = Js.Jv;
        spmd
        Jv = sumJacs(tJv,sJv);
        end
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

spmd
if outputLambda
    [~,g,~,lambdaX,lambdaV] = runSimulateSystemZ(u,x,v,ss,simVars,Jx,Jv);
else
    [~,g] = runSimulateSystemZ(u,x,v,ss,simVars,Jx,Jv);
end

gradU = catAndSum(g);
gradU = gop(@plus,gradU);
end
gradU = gradU{1};

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

JacTarW = cell(size(ss));

JacTarW = cellfun(@(JW,J)subsasgn(JW,struct('type','.','subs','Jx'),J),JacTarW,Jx,'UniformOutput',false);
JacTarW = cellfun(@(JW,J)subsasgn(JW,struct('type','.','subs','Jv'),J),JacTarW,Jv,'UniformOutput',false);

if  nargout >3;

    [f,g,usliced,lambdaX,lambdaV] =...
    cellfun(@(...
    xr,vr,ssr,simVarsr,JacTarWr)...
        simulateSystemZ(u,xr,vr,ssr,[],'simVars',simVarsr,'JacTar',JacTarWr,'withAlgs',true,'printCounter',false),...
    x ,v ,ss ,simVars ,JacTarW ,'UniformOutput',false);
else
    
    [f,g,usliced] =...
        cellfun(@(...
        xr,vr,ssr,simVarsr,JacTarWr)...
        simulateSystemZ(u,xr,vr,ssr,[],'simVars',simVarsr,'JacTar',JacTarWr,'withAlgs',true,'printCounter',false),...
        x ,v ,ss ,simVars ,JacTarW ,'UniformOutput',false);
    
    lambdaX = [];
    lambdaV = [];

end

end
