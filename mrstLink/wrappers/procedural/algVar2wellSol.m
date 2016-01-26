function [ wellSol,JacW,nWV] = algVar2wellSol( v,wellSol,varargin)
%
% write the algebraic variables as a wellSol
%

%  v may contain variables not related to the wells at the tail!
%

opt = struct('vScale',[],...
    'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0),...
    'partials',false);% default OW
opt = merge_options(opt, varargin{:});
comp = opt.activeComponents;


nw = numel(wellSol);


if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    
    
    nWV = nw*3;  %qw qo bhp
    
    if ~isempty(opt.vScale)
        v = v(1:nWV).*opt.vScale(1:nWV);
    end
    
    iqOs = nw;
    ibhp = nw + iqOs;
    
    for w = 1:nw
        wellSol(w).qWs = v(w);
        wellSol(w).qOs = v(iqOs+w);
        wellSol(w).bhp = v(ibhp+w);
    end
    
elseif comp.gas && comp.oil && comp.water
    
    nWV = nw*4; %qw qo qg bhp
    
    if ~isempty(opt.vScale)
        v = v(1:nWV).*opt.vScale(1:nWV);
    end
    
    iqOs = nw;
    iqGs = nw + iqOs;
    ibhp = nw + iqGs;
    
    for w = 1:nw
        wellSol(w).qWs = v(w);
        wellSol(w).qOs = v(iqOs+w);
        wellSol(w).qGs = v(iqGs+w);
        wellSol(w).bhp = v(ibhp+w);
    end
    
else
    error('Not implemented for current activeComponents');
end


if opt.partials
    if ~isempty(opt.vScale)
        JacW = bsxfun(@times,speye(nWV),opt.vScale(1:nWV)');
    else
        JacW = speye(nWV);
    end
else
    JacW = [];
end


end