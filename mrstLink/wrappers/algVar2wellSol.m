function [ wellSol ] = algVar2wellSol( v,wellSol,varargin)
%
% write the algebraic variables as a wellSol
%


opt = struct('vScale',[]);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.vScale)
    v = v.*opt.vScale;
end

nw = numel(wellSol);


iqOs = nw;
ibhp = nw + iqOs;

for w = 1:nw
    wellSol(w).qWs = v(w);
    wellSol(w).qOs = v(iqOs+w);
    wellSol(w).bhp = v(ibhp+w);
end



end