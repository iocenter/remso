function [ wellSol ] = algVar2wellSol( v,wellSol,varargin)
%
% write the algebraic variables as a wellSol
%
%  v may contain variables not related to the wells at the tail!
%

opt = struct('vScale',[]);
opt = merge_options(opt, varargin{:});


nw = numel(wellSol);
nv = nw*3;


if ~isempty(opt.vScale)
    v = v(1:nv).*opt.vScale(1:nv);
end



iqOs = nw;
ibhp = nw + iqOs;

for w = 1:nw
    wellSol(w).qWs = v(w);
    wellSol(w).qOs = v(iqOs+w);
    wellSol(w).bhp = v(ibhp+w);
end



end