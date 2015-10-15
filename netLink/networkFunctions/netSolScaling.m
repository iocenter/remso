function [netSolScale] = netSolScaling(netSol,varargin)
%
% fill a netSol mock object with scaling information according to the
% options
%

opt = struct('bhp',barsa,'qWs',meter^3/day,'qOs',meter^3/day,'qGs',meter^3/day,...
    'activeComponents',struct('oil',1,'water',1,'gas',1));
opt = merge_options(opt, varargin{:});

% comp = opt.activeComponents;
% 
% nw = numel(wellSol);
% netSolScale = netSol;
% 
% 
% if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)    
%     for w= 1:nw
%         netSolScale(w).bhp = opt.bhp;
%         netSolScale(w).qWs = opt.qWs;
%         netSolScale(w).qOs = opt.qOs;
%     end
% elseif comp.gas && comp.oil && comp.water
%     
%     for w= 1:nw
%         netSolScale(w).bhp = opt.bhp;
%         netSolScale(w).qWs = opt.qWs;
%         netSolScale(w).qOs = opt.qOs;
%         netSolScale(w).qGs = opt.qGs;
%     end
%     
%     
% else
%     error('Not implemented for current activeComponents');
% end
    netSolScale = netSol;

end