function [ wellSol ] = algVar2wellSol( v,wellSol,varargin)
%
% write the algebraic variables as a wellSol
%


opt = struct('doScale',false,'vScale',[]);
opt = merge_options(opt, varargin{:});

if opt.doScale
    if isempty(opt.vScale)
        [ wellSolScale ] = wellSolScaling( wellSol );
        [ vScale ] = schedule2Controls(wellSolScale);
        v = v.*vScale;
    else
        v = v.*opt.vScale;
    end
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