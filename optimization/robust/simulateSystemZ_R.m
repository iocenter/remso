function varargout= simulateSystemZ_R(u,x,v,ss,JacTar,varargin)
%
%  Run an adjoint simulation on the Z function in order to compute a
%  gradient with respect to target
%
%
opt = struct('simVars',[],'eta',0.9);
opt = merge_options(opt, varargin{:});

%% Process inputs & prepare outputs

simVars = opt.simVars;
if isempty(simVars)
    simVars = cell(size(ss));
end

if isfield(JacTar,'Js')
    lambdaS = JacTar.Js;
    
    [~,Js] = realization2s(x,u,v,ss,'partials',true,'eta',opt.eta,'leftSeed',lambdaS);
    
    if isfield(JacTar,'Jx')
        JacTar.Jx = sumJacs(JacTar.Jx,Js.Jx);
    else
        JacTar.Jx = Js.Jx;
    end
    if isfield(JacTar,'Jv')
        JacTar.Jv = sumJacs(JacTar.Jv,Js.Jv);
    else
        JacTar.Jv = Js.Jv;
    end
    if isfield(JacTar,'Ju')
        JacTar.Ju = sumJacsS(JacTar.Ju,Js.Ju);
    else
        JacTar.Ju = Js.Ju;
    end
else
    lambdaS = [];
end


JacTarW = cell(size(ss));

JacTarW = cellfun(@(JW,J)subsasgn(JW,struct('type','.','subs','Jx'),J),JacTarW,JacTar.Jx','UniformOutput',false);
JacTarW = cellfun(@(JW,J)subsasgn(JW,struct('type','.','subs','Jv'),J),JacTarW,JacTar.Jv','UniformOutput',false);
% obs -> u is not necessary


[f,g,simVars,usliced,lambdaX,lambdaV] =cellfun(@(xr,vr,ssr,simVarsr,JacTarWr)simulateSystemZ(u,xr,vr,ssr,[],'simVars',simVarsr,'JacTar',JacTarWr,'withAlgs',true),x,v,ss,simVars,JacTarW,'UniformOutput',false);



if isfield(JacTar,'Ju')
    gradU = catAndSum([g;{JacTar.Ju}]);
else
	gradU = catAndSum(g);
end


varargout{1} = gradU;
varargout{2} = lambdaX';
varargout{3} = lambdaV';
varargout{4} = lambdaS;


end



function js = sumJacs(JT,J)

js = cellfun(@sumJacsS,JT,J,'UniformOutput',false);

end
function js = sumJacsS(JT,J)

js = cellfun(@plus,JT,J,'UniformOutput',false);

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
