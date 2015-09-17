

function [netSolScale] = netSOlScaling(netSol,varargin)
%% 
% fill a netSol mock object with scaling information according to the options
%
%opt = struct('bhp',barsa,'qWs',meter^3/day,'qOs',meter^3/day);  %TODO: OW
%opt = merge_options(opt, varargin{:});

%nw = numel(wellSol);

netSolScale = 5*barsa;

%netSolScale = wellSol;



%for w= 1:nw
%    wellSolScale(w).bhp = opt.bhp;
%    wellSolScale(w).qWs = opt.qWs;
%    wellSolScale(w).qOs = opt.qOs;
%end




end