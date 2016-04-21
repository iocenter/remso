function [wellSolScale] = wellSolScaling(wellSol,varargin)
%
% fill a wellSol mock object with scaling information according to the
% options
%

opt = struct('bhp',barsa,'qWs',meter^3/day,'qOs',meter^3/day,'qGs',meter^3/day, 'freq', 1, ...
    'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0));
opt = merge_options(opt, varargin{:});

comp = opt.activeComponents;

nw = numel(wellSol);
wellSolScale = wellSol;


if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    
    for w= 1:nw
        wellSolScale(w).bhp = opt.bhp;
        wellSolScale(w).qWs = opt.qWs;
        wellSolScale(w).qOs = opt.qOs;
    end
    
    
elseif comp.gas && comp.oil && comp.water
    
    for w= 1:nw
        wellSolScale(w).bhp = opt.bhp;
        wellSolScale(w).qWs = opt.qWs;
        wellSolScale(w).qOs = opt.qOs;
        wellSolScale(w).qGs = opt.qGs;
    end
    
    
else
    error('Not implemented for current activeComponents');
end





end