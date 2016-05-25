function [netSolScale] = netSolScaling(netSol,varargin)
%
% fill a wellSol mock object with scaling information according to the
% options
%

opt = struct('pressure', barsa, 'freq', 1, 'flow', meter^3/day);    
opt = merge_options(opt, varargin{:});

equip = getEdge(netSol, netSol.Eeqp);
nN = numel(equip);

if ~isempty(opt.flow) && ~isempty(opt.freq)
    netSolScale = [repmat(opt.flow, nN,1);
                   repmat(opt.flow, nN,1);
                   repmat(opt.freq, nN, 1)
                   ];
elseif ~isempty(opt.pressure)
    netSolScale = [repmat(opt.pressure, nN, 1)];
else
     error('Network constraints not implemented equipment which are neither chokes nor ESPs');
end

end