function [ wellSol,netSol,JacW,JacN ] = algVar2mrstAlg( v,wellSol,netSol,varargin)
%
% write the algebraic variables as a wellSol
%


opt = struct('vScale',[],...
    'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0),...
    'partials',false);% default OW
opt = merge_options(opt, varargin{:});

comp = opt.activeComponents;



[ wellSol,JacW,nWV ] = algVar2wellSol( v,wellSol,'vScale',opt.vScale,...
    'activeComponents',comp,...
    'partials',opt.partials);
nNV = numel(v)-nWV;

%% you may create a similar function algVar2netSol ! for this part

if ~isempty(opt.vScale)
    vN = v(nWV+1:end).*opt.vScale(nWV+1:end);
end

netSol = v(nWV+1:end);

if opt.partials
    if ~isempty(opt.vScale)
        JacN = bsxfun(@times,speye(nNV),opt.vScale(nWV+1:end)');
    else
        JacN = speye(nNV);
    end
else
    JacN = [];    
end



end