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

if isfield(wellSol,'bhp') % TODO: Confirm with MRST Developers bhp or pressure?
    for w = 1:nw
        wellSol(w).qWs = v(w);
        wellSol(w).qOs = v(iqOs+w);
        wellSol(w).bhp = v(ibhp+w);
    end
else
    for w = 1:nw
        wellSol(w).qWs = v(w);
        wellSol(w).qOs = v(iqOs+w);
        wellSol(w).pressure = v(ibhp+w);
    end
end


end