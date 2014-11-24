function [ wellSol,netSol,JacW,JacN ] = algVar2mrstAlg( v,wellSol,netSol,varargin)
%
% write the algebraic variables as a wellSol
%


opt = struct('vScale',[],'partials',false);
opt = merge_options(opt, varargin{:});


if ~isempty(opt.vScale)
    v = v.*[opt.vScale];
end

nw = numel(wellSol);


nWV = nw*3;
nNV = numel(v)-nWV;


iqOs = nw;
ibhp = nw + iqOs;

for w = 1:nw
    wellSol(w).qWs = v(w);
    wellSol(w).qOs = v(iqOs+w);
    wellSol(w).bhp = v(ibhp+w);
end

netSol = v(nWV+1:end);


if opt.partials
    if ~isempty(opt.vScale)
        JacW = bsxfun(@times,speye(nWV),opt.vScale(1:nWV)');
        JacN = bsxfun(@times,speye(nNV),opt.vScale(nWV+1:end)');
    else
        JacW = speye(nWV);
        JacN = speye(nNV);
    end
else
    JacW = [];
    JacN = [];    
end



end