function [ v ] = wellSol2algVar( wellSol,varargin )
%
%  Extract the algebraic variables from the well.
%

opt = struct('vScale',[],...
    'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0));

opt = merge_options(opt, varargin{:});


comp = opt.activeComponents;

if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    
    v = [vertcat(wellSol.qWs);
        vertcat(wellSol.qOs);
        vertcat(wellSol.bhp)];
    
elseif comp.gas && comp.oil && comp.water
    
    v = [vertcat(wellSol.qWs);
        vertcat(wellSol.qOs);
        vertcat(wellSol.qGs);
        vertcat(wellSol.bhp)];
    
else
    error('Not implemented for current activeComponents');
end

nv = numel(v);

if ~isempty(opt.vScale)
    v = v./opt.vScale(1:nv);
end



end

