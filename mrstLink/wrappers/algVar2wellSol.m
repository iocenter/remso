function [ wellSol ] = algVar2wellSol( v,wellSol,varargin)
%
% write the algebraic variables as a wellSol
%


opt = struct('vScale',[],...
    'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0));% default OW
opt = merge_options(opt, varargin{:});

comp = opt.activeComponents;

if ~isempty(opt.vScale)
    v = v.*opt.vScale;
end

nw = numel(wellSol);


if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    
    iqOs = nw;
    ibhp = nw + iqOs;
    
    for w = 1:nw
        wellSol(w).qWs = v(w);
        wellSol(w).qOs = v(iqOs+w);
        wellSol(w).bhp = v(ibhp+w);
    end
    
elseif comp.gas && comp.oil && comp.water
    
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



end